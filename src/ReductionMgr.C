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

 Assumes that *only* one thread will require() a specific sequence's data.
 ***************************************************************************/

static char ident[] = "@(#)$Header: /home/cvs/namd/cvsroot/namd2/src/ReductionMgr.C,v 1.1024 1998/03/03 23:05:26 brunner Exp $";

#include <stdlib.h>
#include <stdio.h>

#include "charm++.h"

#include "InfoStream.h"
#include "PatchMap.h"	// for patchMap

#include "Node.h"
#include "SimParameters.h"


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

/*******************************************
 ReductionMgr::ReductionMgr(): init object
 *******************************************/
ReductionMgr::ReductionMgr(InitMsg *msg)
{
    delete msg;
    if (CpvAccess(ReductionMgr_instance) == 0) {
      CpvAccess(ReductionMgr_instance) = this;
    } else {
      DebugM(1, "ReductionMgr::ReductionMgr() - another instance exists!\n");
    }

    nextSequence = -1; // checked later and changed to firstTimeStep
    data = NULL;
    // data = createdata();

    // fill in the spanning tree fields
    if (CMyPe() == 0) {
      myParent = -1;
    } else {
      myParent = (CMyPe()-1)/MAX_CHILDREN;
    }
    numChildren = 0;
    for(int i=0; i<MAX_CHILDREN; i++) {
      myChildren[i] = CMyPe()*MAX_CHILDREN+i+1;
      if(myChildren[i] < CNumPes())
        numChildren++;
    }

    // initialize data
    for(i=0; i<REDUCTION_MAX_RESERVED; i++)
    {
      numSubscribed[i] = 0;
      maxData[i] = numChildren;
    }

    maxEvents = (numChildren*REDUCTION_MAX_RESERVED); // num events (submit)

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
void ReductionMgr::displayData(ReductionMgrData *current)
{
  if (!current) return;
} /* ReductionMgr::displayData() */

/*******************************************
 ReductionMgr::createdata(): create a blank
 data element.
 *******************************************/
ReductionMgrData *ReductionMgr::createdata()
{
  if ( nextSequence < 0 )
  {
    nextSequence = Node::Object()->simParameters->firstTimestep;
  }

  ReductionMgrData *data;
  data = new ReductionMgrData;
  data->sequenceNum = nextSequence;
  data->next = NULL;
  data->numEvents = 0;
  for(int i=0; i<REDUCTION_MAX_RESERVED; i++)
  {
      data->numData[i] = 0;
      data->tagData[i] = 0;
      data->suspendFlag[i] = 0;
      data->threadNum[i] = 0;
  }
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
 that counts down to zero.  Then it hits 0 it
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

  DebugM(5,"remove()? seq=" << seq << "\n");

  // check if data is removable
  if (currentdata->numEvents < maxEvents)
	return;	// don't delete it!  still required by someone.

  // find data to remove
  for(i=data->sequenceNum; i<seq; i++)
  {
	previousdata = currentdata;
	currentdata = currentdata->next;
  }

  /* don't remove suspended data */
  // for(i=0; i < REDUCTION_MAX_RESERVED; i++)
  // {
  //   if (currentdata->suspendFlag[i]) return;
  // }

  DebugM(5,"remove()! seq=" << seq << "\n");

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
  int seq = msg->seq;
  ReductionMgrData *current=find(seq);
  for(int tag=0;tag<REDUCTION_MAX_RESERVED;tag++) {
    current->tagData[tag] += msg->data[tag];
    current->numData[tag]++;
    current->numEvents++;
  }
  delete msg;

  DebugM(2,"ReductionDataMsg received tag=" << tag
	 << " data=" << current->tagData[tag] << "\n"); 

  // inform object that new data has been found

  DebugM(4,"recv seq=" << current->sequenceNum << " tag=" << tag
	<< " " << current->numData[tag] << "/" << maxData[tag]
	<< " " << current->numEvents << "/" << maxEvents
	<< " data=" << current->tagData[tag]
	<< "\n");

  if (current->numEvents >= maxEvents && !isRoot()) {
    ReductionDataMsg *m 
      = new (MsgIndex(ReductionDataMsg)) ReductionDataMsg;
    m->seq = seq;
    for(tag=0;tag<REDUCTION_MAX_RESERVED;tag++)
      m->data[tag] = current->tagData[tag];
    CSendMsgBranch(ReductionMgr, recvReductionData, ReductionDataMsg, m, thisgroup, myParent);
    gotAllData(current);
  }
  if(isRoot()) {
    for(tag=0;tag<REDUCTION_MAX_RESERVED;tag++) {
      // check if someone is waiting for the data
      // displayData(current,tag);
      if (current->suspendFlag[tag])
      {
          current->suspendFlag[tag] = 0;
          DebugM(5,"Awaken seq=" << current->sequenceNum << " tag=" << tag
	    << " thread=" << current->threadNum[tag]
	    << "\n");
          CthAwaken(current->threadNum[tag]);
      }
      else if (current->numEvents >= maxEvents) {
	    remove(current->sequenceNum);
            break;
      }
    }
  }
} /* ReductionMgr::recvReductionData() */

/*******************************************
 ReductionMgr::find(): find data for reduction.
 If sequence does not exist, then create it.
 *******************************************/
ReductionMgrData *	ReductionMgr::find(int seq)
{
  ReductionMgrData *current=NULL;
  ReductionMgrData *previous=NULL;

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
      DebugM(1,"browsing over seq=" << current->sequenceNum << "\n");
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

  // one less data to send (delete if all done)

  if (current->numEvents >= maxEvents)
	remove(current->sequenceNum);
} /* ReductionMgr::gotAllData() */

/*******************************************
 ReductionMgr::submit(): submit data for reduction.
 more == 1 signals immediate submission of other data
 There should be 1 submit per register.
 *******************************************/
void	ReductionMgr::submit(int seq, ReductionTag tag, BigReal data)
{
  ReductionMgrData *current=find(seq);

  // add to tag
  current->tagData[tag] += data;
  current->numData[tag]++;	/* expect 1 less */
  current->numEvents++;		// got an event (submit)

  DebugM(4,"Submit seq=" << seq
	<< " tag=" << tag
	<< " " << current->numData[tag] << "/" << maxData[tag]
	<< " " << current->numEvents << "/" << maxEvents
	<< " data=" << data
	<< "\n");

  if (current->numEvents >= maxEvents && !isRoot()) {
    ReductionDataMsg *m 
      = new (MsgIndex(ReductionDataMsg)) ReductionDataMsg;
    m->seq = seq;
    for(int i=0;i<REDUCTION_MAX_RESERVED;i++)
      m->data[i] = current->tagData[i];
    CSendMsgBranch(ReductionMgr, recvReductionData, ReductionDataMsg, m, thisgroup, myParent);
    gotAllData(current);
  }
  if (isRoot() && current->numData[tag] == maxData[tag])
  {
	// displayData(current,tag);
	// check if Node 0 (the collector) is suspended
	if (current->suspendFlag[tag])
	{
	  current->suspendFlag[tag] = 0;
	  CthAwaken(current->threadNum[tag]);
	  return;
	}
	else gotAllData(current);
  }

} /* ReductionMgr::submit() */

/*******************************************
 ReductionMgr::submit(): submit data for reduction.
 more == 1 signals immediate submission of other data
 There should be 1 submit per register.
 This function is used when there is NO data to submit.
 *******************************************/
void	ReductionMgr::submit(int seq, ReductionTag tag)
{
  submit(seq,tag,0.0);
} /* ReductionMgr::submit() */

/*******************************************
 ReductionMgr::require(): get data.
 Note: this suspends until this data is ready
 and should be called only from Sequencer thread
 *******************************************/
void	ReductionMgr::require(int seq, ReductionTag tag, BigReal &data)
{
  // 1. find the correct sequence
  ReductionMgrData *current = find(seq);

  // 2. check if all the data is present
  while(current->numData[tag] < maxData[tag])
  {
    // suspend thread until numData is 0...
    current->suspendFlag[tag] = 1;
    current->threadNum[tag] = CthSelf();
    DebugM(5,"Suspend seq=" << seq << " tag=" << tag
	<< " thread=" << current->threadNum[tag]
	<< " " << current->numData[tag] << "/" << maxData[tag]
	<< " " << current->numEvents << "/" << maxEvents
	<< "\n");
    while(current->suspendFlag[tag] == 1)
    {
	current->numEvents--;	// one less event.  Must wait for awaken
	CthSuspend();
	current->numEvents++;	// one more event.  Got awaken
    }
    // ...then return value
    DebugM(5,"unSuspend seq=" << seq << " tag=" << tag
	<< " thread=" << current->threadNum[tag]
	<< "\n");
  }

  // 3. use the data
  data = current->tagData[tag];
  current->numEvents++;		// got an event (require)

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
  numSubscribed[tag]++;
  maxEvents++;
} /* ReductionMgr::subscribe() */

/*******************************************
 ReductionMgr::unsubscribe(): Allow a process
 to no-longer require the data.
 *******************************************/
void	ReductionMgr::unsubscribe(ReductionTag tag)
{
  numSubscribed[tag]--;
  maxEvents--;
} /* ReductionMgr::unsubscribe() */

/*******************************************
 *******************************************/

#include "ReductionMgr.bot.h"
// nothing should be placed below here


/***************************************************************************
 * RCS INFORMATION:
 *
 *	$RCSfile $
 *	$Author $	$Locker:  $		$State: Exp $
 *	$Revision: 1.1024 $	$Date: 1998/03/03 23:05:26 $
 *
 ***************************************************************************
 * REVISION HISTORY:
 *
 * $Log: ReductionMgr.C,v $
 * Revision 1.1024  1998/03/03 23:05:26  brunner
 * Changed include files for new simplified Charm++ include file structure.
 *
 * Revision 1.1023  1998/02/26 01:51:25  milind
 * Fixed bugs in CollectionMaster and ReductionManager that were causing
 * crash on Origin2000.
 *
 * Revision 1.1022  1998/02/10 23:30:32  milind
 * Fixed to reflect the current changes to Charm++ translator.
 *
 * Revision 1.1021  1997/11/07 20:17:48  milind
 * Made NAMD to run on shared memory machines.
 *
 * Revision 1.1020  1997/09/28 10:19:08  milind
 * Fixed priorities, ReductionMgr etc.
 *
 * Revision 1.1019  1997/08/26 16:26:17  jim
 * Revamped prioritites for petter performance and easier changes.
 *
 * Revision 1.1018  1997/08/22 20:12:04  milind
 * Turned on Priorities.
 *
 * Revision 1.1017  1997/07/08 15:48:12  milind
 * Made namd2 to work with Origin2000: Again...
 *
 * Revision 1.1016  1997/04/08 07:08:59  ari
 * Modification for dynamic loadbalancing - moving computes
 * Still bug in new computes or usage of proxies/homepatches.
 * Works if ldbStrategy is none as before.
 *
 * Revision 1.1015  1997/04/07 14:54:35  nealk
 * Changed fclose() to Fclose() (found in common.[Ch]) to use with popen().
 * Also corrected compilation warnings in Set.[Ch].
 *
 * Revision 1.1014  1997/04/06 22:45:13  ari
 * Add priorities to messages.  Mods to help proxies without computes.
 * Added quick enhancement to end of list insertion of ResizeArray(s)
 *
 * Revision 1.1013  1997/04/04 17:31:42  brunner
 * New charm fixes for CommunicateConverse, and LdbCoordinator data file
 * output, required proxies, and idle time.
 *
 * Revision 1.1012  1997/04/02 21:29:41  jim
 * Fixed bad assumption about first sequence number being 0.
 *
 * Revision 1.1011  1997/03/19 11:54:53  ari
 * Add Broadcast mechanism.
 * Fixed RCS Log entries on files that did not have Log entries.
 * Added some register variables to Molecule and ComputeNonbondedExcl.C
 *
 ***************************************************************************/
