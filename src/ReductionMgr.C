/***************************************************************************/
/*                                                                         */
/*              (C) Copyright 1996 The Board of Trustees of the            */
/*                          University of Illinois                         */
/*                           All Rights Reserved                           */
/*                                                                         */
/***************************************************************************/

/***************************************************************************
 * DESCRIPTION:
 *
 ***************************************************************************/

/***************************************************************************
 The order of execution is expected to be:
	    0. instantiate object
   ------------------ (processing barrier)
   (mode 0) 1. register() and subscribe()
   ------------------
   (mode 1) 2. submit() and request()
   ------------------
   (mode 2) 3. unregister() and unsubscribe()
   ------------------ (processing barrier)
            4. destroy object
 Doing this out-of-order will cause errors.

 Assumes that *only* node 0 will require data.
 Assumes that *only* one thread will require() a specific sequence's data.
 ***************************************************************************/

static char ident[] = "@(#)$Header: /home/cvs/namd/cvsroot/namd2/src/ReductionMgr.C,v 1.781 1997/02/06 15:53:27 ari Exp $";

#include <stdlib.h>
#include <stdio.h>

#include "chare.h"
#include "ckdefs.h"
#include "c++interface.h"

#include "InfoStream.h"
#include "PatchMap.h"	// for patchMap

// determine whether PANIC sequence checking is performed (debugging)
#define PANIC 0
#include "ReductionMgr.top.h"
#include "ReductionMgr.h"

// #define DEBUGM
#define MIN_DEBUG_LEVEL 4
#define STDERR_LEVEL 7
#include "Debug.h"

// *************************** for Charm messages

/** questions:
 1. who is responsible for allocating/deleting the charm data from pack/unpack
 2. how is the thread "awakened"
 3. who activates the pack/unpack functions
**/

// *************************** for object control

ReductionMgr *ReductionMgr::_instance = 0;

/*******************************************
 ReductionMgr::ReductionMgr(): init object
 *******************************************/
ReductionMgr::ReductionMgr(InitMsg *msg)
{
    #if PANIC > 0
    panicMode = 0;
    strcpy(tagString[0],"REDUCTION_ANGLE_ENERGY");
    strcpy(tagString[1],"REDUCTION_BOND_ENERGY");
    strcpy(tagString[2],"REDUCTION_DIHEDRAL_ENERGY");
    strcpy(tagString[3],"REDUCTION_ELECT_ENERGY");
    strcpy(tagString[4],"REDUCTION_IMPROPER_ENERGY");
    strcpy(tagString[5],"REDUCTION_KINETIC_ENERGY");
    strcpy(tagString[6],"REDUCTION_LJ_ENERGY");
    strcpy(tagString[7],"REDUCTION_LONG_RANGE_ENERGY");
    strcpy(tagString[8],"REDUCTION_MAX_RESERVED");
    #endif

    delete msg;
    if (_instance == 0) {
      _instance = this;
    } else {
      DebugM(1, "ReductionMgr::ReductionMgr() - another instance exists!\n");
    }

    nextSequence = 0;
    data = createdata();

    /* node 0 received data from each remove node */
    if (CMyPe() == 0)
    {
      // Problem: with small systems, there may be fewer patches than nodes
      // This problem causes Node 0 to wait for data that never arrives.
      // (Nodes without patches never compute!)

      // node 0 does not submit to node 0 (but every other Pe does)
      int num = CNumPes()-1;
      maxEvents += (num*REDUCTION_MAX_RESERVED); // num events (submit)

      DebugM(1,"Initializing with a minimum of " << num << " data.\n");
      for(int i=0; i<REDUCTION_MAX_RESERVED; i++)
      {
    	numSubscribed[i] = 0;
	maxData[i] = num;
      }
    }
    else
    {
      maxEvents = 0;
      for(int i=0; i<REDUCTION_MAX_RESERVED; i++)
      {
    	numSubscribed[i] = 0;
	maxData[i] = 0;
      }
    }

    DebugM(1,"ReductionMgr() instantiated.\n");
} /* ReductionMgr::ReductionMgr() */

/*******************************************
 ReductionMgr::~ReductionMgr(): free object
 *******************************************/
ReductionMgr::~ReductionMgr()
{
    ReductionMgrData *nextdata;
    while(data != NULL)
    {
      nextdata = data->next;
      delete data;
      data = nextdata;
    }
} /* ReductionMgr::~ReductionMgr() */

/*******************************************
 ReductionMgr::displayData(): dump the data.
 Very useful for debugging.
 *******************************************/
void ReductionMgr::displayData(ReductionMgrData *current, ReductionTag tag)
{
   #if PANIC > 0
   iout << iPE << " seq=" << current->sequenceNum
	<< " " << tagString[tag]
	<< " " << current->tagData[tag] << "\n" << endi;
   #endif
} /* ReductionMgr::displayData() */

/*******************************************
 ReductionMgr::displayData(): dump the data.
 Very useful for debugging.
 *******************************************/
void ReductionMgr::displayData(ReductionMgrData *current)
{
  if (!current) return;
  #if PANIC > 0
  if (current->dataToSend)
	iout << "Unfilled data fields: "<< current->dataToSend <<"\n" << endi;
  for(int tag=0; tag<REDUCTION_MAX_RESERVED; tag++)
	displayData(current,(ReductionTag)tag);
  #endif
} /* ReductionMgr::displayData() */

/*******************************************
 ReductionMgr::createdata(): create a blank
 data element.
 *******************************************/
ReductionMgrData *ReductionMgr::createdata()
{
  ReductionMgrData *data;
  data = new ReductionMgrData;
  data->sequenceNum = nextSequence;
  data->next = NULL;
  data->numEvents = 0;
  #if PANIC > 0
  data->dataToSend = REDUCTION_MAX_RESERVED;
  #endif

  for(int i=0; i<REDUCTION_MAX_RESERVED; i++)
  {
      data->numData[i] = 0;
      data->tagData[i] = 0;
      data->suspendFlag[i] = 0;
      data->threadNum[i] = 0;
  }
  data->eventCounter = 0;

  DebugM(4," createdata(" << nextSequence << ")\n");
  nextSequence++;
  return(data);
} /* ReductionMgr::createdata() */

/*******************************************
  ReductionMgr::Register(): increase counter
  to a reduction tag.
  Only registered objects may deposit data.

  ASSUMPTION: this function will only be called
  before data has been depositied.
  (un)register to submit data for reduction
  may cause an error if reductions are active
 *******************************************/
void	ReductionMgr::Register(ReductionTag tag)
{
  #if PANIC > 0
  if (panicMode != 0)
  {
    iout << iERRORF << "Panic due to wrong mode: " << panicMode << "\n" << endi;
    NAMD_die("Panic due to wrong mode");
  }
  panicMode = 0;
  #endif

  maxData[tag]++;
  maxEvents++;	// expect and event (submit)
  DebugM(1,"Register tag=" << tag << " maxData="<< maxData[tag] <<"\n");
} /* ReductionMgr::Register() */

/*******************************************
 ReductionMgr::unRegister(): 
 ASSUMPTION: this function will only be called
 after data has been depositied.
 (un)register to submit data for reduction
 may cause an error if reductions are active
 *******************************************/
void	ReductionMgr::unRegister(ReductionTag tag)
{
  #if PANIC > 0
  if (panicMode < 1)
  {
    iout << iERRORF << "Panic due to wrong mode: " << panicMode << "\n" << endi;
    NAMD_die("Panic due to wrong mode");
  }
  panicMode = 2;
  #endif

  maxData[tag]--;
  maxEvents--;	// expect 1 less event
  DebugM(1,"unRegister tag=" << tag << " maxData=" << maxData[tag] << "\n");
} /* ReductionMgr::unRegister() */

/*******************************************
 ReductionMgr::remove(): remove a sequence
 of counters.
 ASSUMPTION: this function will only be called
 after data has been depositied.
 (un)register to submit data for reduction
 may cause an error if reductions are active

 Note: there should be a general event flag
 that counts down to zero.  When it hits 0 it
 should remove the sequence.
   numEvents = numRequire+numSubmit
 This will cause a minor performance improval.
 Update: Now we count events!
 *******************************************/
void	ReductionMgr::remove(int seq)
{
  if (!data) return;
  ReductionMgrData *currentdata = data;
  ReductionMgrData *previousdata = NULL;
  int i;	// loop variable

  DebugM(4,"remove()? seq=" << seq << "\n");

  // check if data is removable
  if (currentdata->numEvents < maxEvents)
	return;	// don't delete it!  still required by someone.

  // find data to remove
  for(i=data->sequenceNum; i<seq; i++)
  {
	#if PANIC > 0
	if (!currentdata)
	  {
	  iout << iERRORF << "Oh no, Mr. Bill!  Cannot find " << seq
	       << " to remove!\n" << endi;
	  NAMD_die("Cannot remove sequence.");
	  }
	#endif
	previousdata = currentdata;
	currentdata = currentdata->next;
  }

  #if PANIC > 0
  if (currentdata->sequenceNum != seq)
  {
    iout << iERRORF << "Yikes! Missed sequence " << seq << " and found "
	 << currentdata->sequenceNum << " instead!\n" << endi;
    NAMD_die("Remove without register");
  }
  #endif

  DebugM(4,"remove()! seq=" << seq << "\n");

  // delete data
  if (!previousdata)
    {
      // head of queue
      data = currentdata->next;
    }
  else
    {
      // body of queue
      previousdata = currentdata->next;
    }
  delete currentdata;
} /* ReductionMgr::remove() */

/*******************************************
 ReductionMgr::recvReductionData(): receive
 and include some data from a BOC.
 *******************************************/
void	ReductionMgr::recvReductionData	(ReductionDataMsg *msg)
{
  ReductionMgrData *current=find(msg->seq);
  ReductionTag tag = msg->tag;
  current->tagData[tag] += msg->data;
  delete msg;

  #if PANIC > 0
  DebugM(2,"ReductionDataMsg received tag=" << tagString[tag]
	 << " data=" << current->tagData[tag] << "\n"); 
  #else
  DebugM(2,"ReductionDataMsg received tag=" << tag
	 << " data=" << current->tagData[tag] << "\n"); 
  #endif

  // inform object that new data has been found
  current->numData[tag]++;
  current->numEvents++;

  DebugM(4,"recv seq=" << current->sequenceNum << " tag=" << tag
	<< " data=" << current->tagData[tag]
	<< " #" << current->numEvents << "/" << maxEvents << "\n");

  #if PANIC > 0
  if (current->numData[tag] > maxData[tag])
  {
    iout << iPE << " " << iERRORF << "Too many receives: seq="
	 << current->sequenceNum
	 << " tag=" << tagString[tag]
	 << " sum+data=" << current->tagData[tag]
	 << "\n" << endi;
  }
  else if (current->numData[tag] == maxData[tag])
  {
    DebugM(3, " Reduction data for seq=" << current->sequenceNum
	 << " tag=" << tagString[tag]
	 << " is " << current->tagData[tag] << "\n");
  }
  #endif

  // check if someone is waiting for the data
  if (current->numData[tag] == maxData[tag])
  {
    // displayData(current,tag);
    if (current->suspendFlag[tag])
    {
      current->suspendFlag[tag] = 0;
      DebugM(5,"Awaken seq=" << current->sequenceNum << " tag=" << tag
	<< " thread=" << current->threadNum[tag]
	<< "\n");
      CthAwaken(current->threadNum[tag]);
    }
    gotAllData(current);
  }

  if (current->numEvents >= maxEvents)
	remove(current->sequenceNum);
} /* ReductionMgr::recvReductionData() */

/*******************************************
 ReductionMgr::find(): find data for reduction.
 If sequence does not exist, then create it.
 *******************************************/
ReductionMgrData *	ReductionMgr::find(int seq)
{
  ReductionMgrData *current=NULL;
  ReductionMgrData *previous=NULL;

  #if PANIC > 0
  if (data && (seq < data->sequenceNum))
  {
    iout << iPE << iERRORF << "Searching for data that has been removed. "
	 << "want " << seq << " but queue starts at " << data->sequenceNum
	 << "\n" << endi;
    NAMD_die("Searching for data that has been removed");
  }
  #endif

  // check for an empty queue
  if (!data)
  {
    data = createdata();
  }
  current = data;

  // find the sequence
  while(current && (current->sequenceNum < seq))
    {
      previous = current;
      DebugM(2,"browsing over seq=" << current->sequenceNum << "\n");
      current = current->next;
    }

  // check if sequence needs to be added
  if (!current)
  {
    while(previous->sequenceNum < seq)
    {
      previous->next = createdata();
      previous->next->next = NULL;
      previous = previous->next;
    }
    current = previous;
  }

  // current now contains the sequence
  #if PANIC > 0
  if (current->sequenceNum != seq)
  {
	iout << iPE << " " << iERRORF
	     << "Current queue has wrong sequence number: seq=" << seq
	     << " current->sequenceNum=" << current->sequenceNum
	     << " Qstart=" << data->sequenceNum
	     << "\n" << endi;
	NAMD_die("Panic due to wrong sequence");
  }
  #endif

  return(current);
} /* ReductionMgr::find() */

/*******************************************
 ReductionMgr::gotAllData(): things to do when
 all data in a sequence is received.
 Currently does:
   1. prints to terminal
   2. frees memory
 If not all data has been received, then returns quietly
 *******************************************/
void	ReductionMgr::gotAllData(ReductionMgrData *current)
{
  DebugM(2,"All data collected for seq=" << current->sequenceNum << "\n");

  #if PANIC > 0
  // one less data to send (delete if all done)
  #if PANIC > 0
  current->dataToSend--;
  if (current->dataToSend < 0)
  {
    iout << iERRORF << iPE << " gotAllData(): received too much."
	 << " seq=" << current->sequenceNum
	 << "\n";
  }

  // display all the data
  if (current->dataToSend == 0)
  {
    displayData(current);
  }

  // remove when done
  if (current->dataToSend <= 0)
  {
      if (CMyPe() == 0)
      {
	DebugM(4,"Node 0 reduction done.\n");
      }
  }
  #endif

  if (current->numEvents >= maxEvents)
	remove(current->sequenceNum);
} /* ReductionMgr::gotAllData() */

/*******************************************
 ReductionMgr::submit(): submit data for reduction.
 more == 1 signals immediate submission of other data
 There should be 1 submit per register.
 *******************************************/
void	ReductionMgr::submit(int seq, ReductionTag tag, BigReal data, int more)
{
  #if PANIC > 0
  if (panicMode > 1)
  {
    iout << iERRORF << "Panic due to wrong mode: " << panicMode << "\n" << endi;
    NAMD_die("Panic due to wrong mode");
  }
  panicMode = 1;
  #endif

  ReductionMgrData *current=find(seq);

  // add to tag
  current->tagData[tag] += data;
  current->numData[tag]++;	/* expect 1 less */
  current->numEvents++;		// got an event (submit)

  DebugM(1,"Submit seq=" << seq << " tag=" << tag << " data=" << data
	<< " #" << current->numData[tag] << "/" << maxData[tag]
	<< " #" << current->numEvents << "/" << maxEvents << "\n");

  #if PANIC > 0
  if (current->numData[tag] > maxData[tag])
  {
    iout << iPE << " " << iERRORF << "Too many submits: seq=" << seq
	 << " tag=" << tagString[tag]
	 << " data=" << data
	 << "\n" << endi;
  }
  else if (current->numData[tag] == maxData[tag])
  {
    DebugM(3, " Reduction data for seq=" << seq << " tag=" << tagString[tag]
	 << " is " << current->tagData[tag] << "\n");
  }
  #endif

  if (current->numData[tag] == maxData[tag])
  {
    // send the data[tag] to main collector
    // don't send the data if this object IS the collector
    if (CMyPe() != 0)
    {
      ReductionDataMsg *m = new (MsgIndex(ReductionDataMsg)) ReductionDataMsg;
      m->seq = seq;
      m->tag = tag;
      m->data = current->tagData[tag];
      DebugM(4,"Sending seq=" << seq << " tag=" << tag
     		 << " data=" << m->data << "\n");
      CSendMsgBranch(ReductionMgr, recvReductionData, m, thisgroup, 0);
      DebugM(3,"Sent seq=" << seq << " tag=" << tag
     		 << " data=" << current->tagData[tag] << "\n");
    }
    else
    {
      if (current->numData[tag] == maxData[tag])
      {
	// displayData(current,tag);
	// check if Node 0 (the collector) is suspended
	if (current->suspendFlag[tag])
	{
	  current->suspendFlag[tag] = 0;
	  CthAwaken(current->threadNum[tag]);
	}
      }
    }

    // all done here!
    gotAllData(current);
  }
} /* ReductionMgr::submit() */

/*******************************************
 ReductionMgr::submit(): pass on submitting data
 *******************************************/
void	ReductionMgr::submit(int seq, ReductionTag tag, int more)
{
  #if PANIC > 0
  if (panicMode > 1)
  {
    iout << iERRORF << "Panic due to wrong mode: " << panicMode << "\n" << endi;
    NAMD_die("Panic due to wrong mode");
  }
  panicMode = 1;
  #endif
  DebugM(2,"Other submit tag=" << tag << "\n");
} /* ReductionMgr::submit() */

/*******************************************
 ReductionMgr::require(): get data.
 Note: this suspends until this data is ready
 and should be called only from Sequencer thread
 *******************************************/
void	ReductionMgr::require(int seq, ReductionTag tag, BigReal &data)
{
  #if PANIC > 0
  if (panicMode > 1)
  {
    iout << iERRORF << "Panic due to wrong mode: " << panicMode << "\n" << endi;
    NAMD_die("Panic due to wrong mode");
  }
  panicMode = 1;
  #endif

  // 1. find the correct sequence
  ReductionMgrData *current = find(seq);

  // 2. check if all the data is present
  while(current->numData[tag] < maxData[tag])
  {
    #if PANIC > 0
    if (current->suspendFlag[tag])
    {
      iout << iERRORF << iPE << " two threads are waiting for the same data! "
	   << current->threadNum << " is suspended and " << CthSelf()
	   << " wants to suspend.\n" << endi;
      NAMD_die("Two threads are waiting for the same data");
    }
    #endif

    // suspend thread until numData is 0...
    current->suspendFlag[tag] = 1;
    current->threadNum[tag] = CthSelf();
    DebugM(5,"Suspend seq=" << seq << " tag=" << tag
	<< " thread=" << current->threadNum[tag]
	<< "\n");
    while(current->suspendFlag[tag] == 1)
	CthSuspend();
    // ...then return value
    if (current->suspendFlag[tag] == 1)
    {
	iout << iERRORF << iPE << " Processes didn't suspend!\n" << endi;
	NAMD_die("ReductionMgr.C: Processes didn't suspend!");
    }
    DebugM(5,"unSuspend seq=" << seq << " tag=" << tag
	<< " thread=" << current->threadNum[tag]
	<< "\n");
  }

  // 3. use the data
  data = current->tagData[tag];
  current->numEvents++;		// got an event (require)
  #if PANIC > 1
  iout << iPE << " Data for sequence=" << seq << " tag=" << tagString[tag] << " is "
       << data << "\n" << endi;
  #endif

  if (current->numEvents == maxEvents)
	remove(seq);    // free it.
} /* ReductionMgr::require() */

/*******************************************
 ReductionMgr::subscribe(): Allow a process to
 require the data.
 A process must subscribe if it requires data.
 *******************************************/
void	ReductionMgr::subscribe(ReductionTag tag)
{
  #if PANIC > 0
  if (panicMode != 0)
  {
    iout << iERRORF << "Panic due to wrong mode: " << panicMode << "\n" << endi;
    NAMD_die("Panic due to wrong mode");
  }
  panicMode = 0;
  #endif

  numSubscribed[tag]++;
  maxEvents++;
} /* ReductionMgr::subscribe() */

/*******************************************
 ReductionMgr::unsubscribe(): Allow a process
 to no-longer require the data.
 *******************************************/
void	ReductionMgr::unsubscribe(ReductionTag tag)
{
  #if PANIC > 0
  panicMode = 2;
  #endif

  numSubscribed[tag]--;
  maxEvents--;
  #if PANIC > 0
  if (numSubscribed[tag] < 0)
  {
    iout << iERRORF << "Too many unsubscribed. " << numSubscribed[tag]
	 << " for tag==" << tagString[tag] << "\n" << endi;
    NAMD_die("Too many unsubscribed.");
  }
  #endif
} /* ReductionMgr::unsubscribe() */

/*******************************************
 *******************************************/

#include "ReductionMgr.bot.h"
// nothing should be placed below here

