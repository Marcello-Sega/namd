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

static char ident[] = "@(#)$Header: /home/cvs/namd/cvsroot/namd2/src/ReductionMgr.C,v 1.12 1997/01/16 21:09:38 nealk Exp $";

#include <stdlib.h>
#include <stdio.h>

#include "chare.h"
#include "ckdefs.h"
#include "c++interface.h"

#include "InfoStream.h"
#include "ReductionMgr.top.h"
#include "ReductionMgr.h"
#define DEBUGM
#define MIN_DEBUG_LEVEL 4
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

    data = NULL;

    /* node 0 received data from each remove node */
    if (CMyPe() == 0)
    {
      for(int i=0; i<REDUCTION_MAX_RESERVED; i++)
      {
    	numSubscribed[i] = 0;
	maxData[i] = CNumPes()-1;
      }
    }
    else
    {
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
      #if PANIC > 0
      displayData(data);
      #endif
      delete data;
      data = nextdata;
    }
} /* ReductionMgr::~ReductionMgr() */

/*******************************************
 ReductionMgr::displayData(): dump the data.
 Very useful for debugging.
 *******************************************/
void ReductionMgr::displayData(ReductionMgrData *current)
{
  if (!current) return;
  if (current->dataToSend)
	iout << "Unfilled data fields: " << current->dataToSend << "\n" << endi;
  for(int tag=0; tag<REDUCTION_MAX_RESERVED; tag++)
      iout << iPE << " " << current->sequenceNum
	   << " " << tagString[tag]
	   << " " << current->tagData[tag] << "\n" << endi;
} /* ReductionMgr::displayData() */

/*******************************************
 ReductionMgr::createdata(): create a blank
 data element.
 *******************************************/
ReductionMgrData *ReductionMgr::createdata(int seq)
{
  ReductionMgrData *data;
  data = new ReductionMgrData;
  data->sequenceNum = seq;
  data->dataToSend = REDUCTION_MAX_RESERVED;
  // data->dataToSend = 1;
  data->next = NULL;
  data->threadNum = 0;
  data->suspendFlag = 0;
  for(int i=0; i<REDUCTION_MAX_RESERVED; i++)
  {
      data->numData[i] = maxData[i];
      data->tagData[i] = 0;
  }
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
  DebugM(1,"Register tag=" << tag << " maxData=" << maxData[tag] << "\n");
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
  DebugM(1,"unRegister tag=" << tag << " maxData=" << maxData[tag] << "\n");
} /* ReductionMgr::unRegister() */

/*******************************************
 ReductionMgr::remove(): remove a sequence
 of counters.
 ASSUMPTION: this function will only be called
 after data has been depositied.
 (un)register to submit data for reduction
 may cause an error if reductions are active
 *******************************************/
void	ReductionMgr::remove(int seq)
{
  if (!data) return;
  ReductionMgrData *currentdata = data;
  ReductionMgrData *previousdata = NULL;

  DebugM(1,"remove() seq=" << seq << "\n");
  for(int i=data->sequenceNum; i<seq; i++)
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

  DebugM(2,"ReductionDataMsg received tag=" << tag << " data="
         << current->tagData[tag] << "\n"); 

  // inform object that new data has been found
  current->numData[tag]--;

  DebugM(1,"recv seq=" << current->sequenceNum << " tag=" << tag
	<< " data=" << current->tagData[tag]
	<< " #" << current->numData[tag] << "/" << maxData[tag] << "\n");

  #if PANIC > 0
  if (current->numData[tag] < 0)
  {
    iout << iPE << " " << iERRORF << "Too many receives: seq="
	 << current->sequenceNum
	 << " tag=" << tag
	 << " sum+data=" << current->tagData[tag]
	 << "\n" << endi;
  }
  else if (current->numData[tag] == 0)
  {
    iout << iPE << " Reduction data for seq=" << current->sequenceNum
	 << " tag=" << tag
	 << " is " << current->tagData[tag] << "\n" << endi;
  }
  #endif

  // check if someone is waiting for the data
  if (current->numData[tag] == 0)
  {
    gotAllData(current);
    if (current->suspendFlag)
    {
      current->suspendFlag = 0;
      CthAwaken(current->threadNum);
    }
  }
} /* ReductionMgr::recvReductionData() */

/*******************************************
 ReductionMgr::find(): find data for reduction.
 If sequence does not exist, then create it.
 *******************************************/
ReductionMgrData *	ReductionMgr::find(int seq)
{
  ReductionMgrData *current=data;
  ReductionMgrData *previous=NULL;

  #if PANIC > 0
  if (current && (seq < current->sequenceNum))
  {
    iout << iERRORF << "Searching for data that has been removed. "
	 << "want " << seq << " but queue starts at " << current->sequenceNum
	 << "\n" << endi;
    NAMD_die("Searching for data that has been removed");
  }
  #endif

  // find the sequence
  while(current && (current->sequenceNum < seq))
  {
    previous = current;
    current = current->next;
  }

  if (!current)
    {
    // no sequences available (create)
    if (!previous)
	{
	data = createdata(seq);
	current = data;
	}
    // it's a middle sequence (insert)
    else
	{
	previous->next = createdata(seq);
	previous->next->next = current;
	current = previous->next;
	}
    }
  // current now contains the sequence
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
  DebugM(2,"All data collected for sequence=" << current->sequenceNum << "\n");

  // one less data to send (delete if all done)
  current->dataToSend--;

  #if PANIC > 0
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
  #endif

  // remove when done
  if (current->dataToSend == 0)
  {
      if (CMyPe() != 0) remove(current->sequenceNum);
      else
      {
	iout << "Node 0 reduction done.\n" << endi;
      }
  }
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
  current->numData[tag]--;	/* expect 1 less */

  DebugM(1,"Submit seq=" << seq << " tag=" << tag << " data=" << data
	<< " #" << current->numData[tag] << "/" << maxData[tag] << "\n");

  #if PANIC > 0
  if (current->numData[tag] < 0)
  {
    iout << iPE << " " << iERRORF << "Too many submits: seq=" << seq
	 << " tag=" << tag
	 << " data=" << data
	 << "\n" << endi;
  }
  else if (current->numData[tag] == 0)
  {
    iout << iPE << " Reduction data for seq=" << seq << " tag=" << tag
	 << " is " << current->tagData[tag] << "\n" << endi;
  }
  #endif

  if (current->numData[tag] == 0)
  {
    // send the data[tag] to main collector
    // don't send the data if this object IS the collector
    if (CMyPe() != 0)
    {
      ReductionDataMsg *m = new (MsgIndex(ReductionDataMsg)) ReductionDataMsg;
      m->seq = seq;
      m->tag = tag;
      m->data = current->tagData[tag];
      CSendMsgBranch(ReductionMgr, recvReductionData, m, thisgroup, 0);
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
  DebugM(1,"Other submit tag=" << tag << "\n");
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
  if (current->numData[tag] > 0)
  {
    #if PANIC > 0
    if (current->suspendFlag)
    {
      iout << iERRORF << iPE << " two threads are waiting for the same data! "
	   << current->threadNum << " is suspended and " << CthSelf()
	   << " wants to suspend.\n" << endi;
      NAMD_die("Two threads are waiting for the same data");
    }
    #endif

    // suspend thread until numData is 0...
    current->suspendFlag = 1;
    current->threadNum = CthSelf();
    DebugM(1,"Thread " << current->threadNum << " waiting for "
			<< current->numData[tag] << " data.\n");
    CthSuspend();
    // ...then return value
  }

  // 3. use the data
  data = current->tagData[tag];
  #if PANIC > 0
  iout << iPE << " Data for sequence=" << seq << " tag=" << tag << " is "
       << data << "\n" << endi;
  #endif

  current->numRequire[tag]--;
  if (current->numRequire[tag] == 0)
	remove(seq);	// free it.
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
  #if PANIC > 0
  if (numSubscribed[tag] < 0)
  {
    iout << iERRORF << "Too many unsubscribed. " << numSubscribed[tag]
	 << " for tag==" << tag << "\n" << endi;
    NAMD_die("Too many unsubscribed.");
  }
  #endif
} /* ReductionMgr::unsubscribe() */

/*******************************************
 *******************************************/

#include "ReductionMgr.bot.h"
// nothing should be placed below here

