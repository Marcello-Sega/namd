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

static char ident[] = "@(#)$Header: /home/cvs/namd/cvsroot/namd2/src/ReductionMgr.C,v 1.5 1996/12/27 22:22:33 nealk Exp $";

#include <stdlib.h>
#include <stdio.h>

#include "chare.h"
#include "ckdefs.h"
#include "c++interface.h"

#include "InfoStream.h"
#include "ReductionMgr.top.h"
#include "ReductionMgr.h"
#define DEBUGM
#include "Debug.h"

// *************************** for Charm messages

/** questions:
 1. who is responsible for allocating/deleting the charm data from pack/unpack
 2. how is the thread "awakened"
 3. who activates the pack/unpack functions
**/

// *************************** for object control

/*******************************************
 ReductionMgr::ReductionMgr(): init object
 *******************************************/
ReductionMgr::ReductionMgr()
{
    #if PANIC > 0
    panicMode = 0;
    #endif
    data = NULL;
    for(int i=0; i<REDUCTION_MAX_RESERVED; i++)
    {
    	numSubscribed[i] = 0;
	maxData[i] = 0;
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
 ReductionMgr::createdata(): create a blank
 data element.
 *******************************************/
ReductionMgrData *ReductionMgr::createdata(int seq)
{
  ReductionMgrData *data;
  data = new ReductionMgrData;
  data->sequenceNum = seq;
  data->dataToSend = REDUCTION_MAX_RESERVED;
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

  // inform object that new data has been found
  current->numData[tag]--;

  // check if someone is waiting for the data
  if ((current->numData[tag] == 0) && (current->suspendFlag))
  {
    current->suspendFlag = 0;
    CthAwaken(current->threadNum);
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

  #if PANIC > 0
  if (current->numData[tag] < 0)
  {
    iout << iERRORF << "Too many submits!\n" << endi;
    NAMD_die("Too many submits");
  }
  else if (current->numData[tag] == 0)
  {
    iout << iPE << " Reduction data for seq=" << seq << " tag=" << tag
	 << " is " << current->tagData[tag] << "\n" << endi;
  }
  #endif

  if (current->numData == 0)
  {
    // all done here!
    // send the data[tag] to main collector
    ReductionDataMsg *m = new ReductionDataMsg;
    m->seq = seq;
    m->tag = tag;
    m->data = current->tagData[tag];
    CSendMsgBranch(ReductionMgr, recvReductionData, m, thisgroup, 0);

    // one less data to send (delete if all done)
    current->dataToSend--;
    if ((CMyPe() != 0) && (current->dataToSend == 0))
	remove(seq);
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
