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

#include <stdlib.h>
#include <stdio.h>

#include "charm++.h"

#include "InfoStream.h"
#include "PatchMap.h"	// for patchMap

#include "Node.h"
#include "SimParameters.h"

#include "ReductionMgr.decl.h"
#include "ReductionMgr.h"

// #define DEBUGM
#define MIN_DEBUG_LEVEL 4
#include "Debug.h"

// Later this can be dynamic
#define REDUCTION_MAX_CHILDREN 4

// Used to register and unregister reductions to downstream nodes
class ReductionRegisterMsg : public CMessage_ReductionRegisterMsg {
public:
  int reductionSetID;
  int sourceNode;
};

// Used to send reduction data to downstream nodes
class ReductionSubmitMsg : public CMessage_ReductionSubmitMsg {
public:
  int reductionSetID;
  int sourceNode;
  int sequenceNumber;
  int dataSize;
  BigReal *data;

  static void *alloc(int msgnum, int size, int *array, int priobits) {
    int totalsize = size + array[0]*sizeof(BigReal);
    ReductionSubmitMsg *newMsg = (ReductionSubmitMsg*)
				CkAllocMsg(msgnum,totalsize,priobits);
    newMsg->data = (BigReal*) ((char*)newMsg + size);
    return (void*)newMsg;
  }

  static void *pack(ReductionSubmitMsg *in) {
    in->data = (BigReal*) ((char*)in->data - (char*)&(in->data));
    return (void*)in;
  }

  static ReductionSubmitMsg *unpack(void *in) {
    ReductionSubmitMsg *me = (ReductionSubmitMsg*)in;
    me->data = (BigReal*) ((char*)&(me->data) + (size_t)(me->data));
    return me;
  }

};

// Queue element which stores data for a particular sequence number
struct ReductionSetData {
  int sequenceNumber;
  int eventsRemaining;  // includes delivery, NOT suspend
  int dataSize;
  BigReal *data;
  ReductionSetData *next;
  ReductionSetData(int seqNum, int events) {
    sequenceNumber = seqNum;
    eventsRemaining = events;
    dataSize = 0;
    data = 0;
    next = 0;
  }
  ~ReductionSetData(void) {
    delete [] data;
  }
  void resize(int size) {
    if ( size > dataSize ) {
      BigReal *oldData = data;
      data = new BigReal[size];
      int i = 0;
      for ( ; i < dataSize; ++i ) { data[i] = oldData[i]; }
      for ( ; i < size; ++i ) { data[i] = 0; }
      dataSize = size;
      delete [] oldData;
    }
  }
};

// Stores the submit queue for a particular set of reductions
struct ReductionSet {
  int reductionSetID;
  int nextSequenceNumber;
  int eventsRegistered;
  ReductionSetData *dataQueue;
  ReductionSetData* getData(int seqNum);
  void delData(int seqNum);
  int requireRegistered;  // is a thread subscribed on this node?
  int threadIsWaiting;  // is there a thread waiting on this?
  int waitingForSequenceNumber;  // sequence number waited for
  CthThread waitingThread;
  ReductionSet(int setID) {
    reductionSetID = setID;
    nextSequenceNumber = 0;
    eventsRegistered = 0;
    dataQueue = 0;
    requireRegistered = 0;
    threadIsWaiting = 0;
  }
  int addToRemoteSequenceNumber[REDUCTION_MAX_CHILDREN];
};

// possibly create and return data for a particular seqNum
ReductionSetData* ReductionSet::getData(int seqNum) {

  ReductionSetData **current = &dataQueue;

  while ( *current ) {
    if ( (*current)->sequenceNumber == seqNum ) return *current;
    current = &((*current)->next);
  }

  nextSequenceNumber++; // should match all clients
  *current = new ReductionSetData(seqNum, eventsRegistered);
  return *current;
}

// possibly delete data for a particular seqNum
void ReductionSet::delData(int seqNum) {

  ReductionSetData **current = &dataQueue;

  while ( *current ) {
    if ( (*current)->sequenceNumber == seqNum ) break;
    current = &((*current)->next);
  }

  if ( ! *current ) { NAMD_die("ReductionSet::delData on missing seqNum"); }

  if ( (*current)->eventsRemaining == 0 ) {
    ReductionSetData *todelete = *current;
    *current = (*current)->next;
    delete todelete;
  }
}

// constructor
ReductionMgr::ReductionMgr() {
    if (CpvAccess(ReductionMgr_instance) == 0) {
      CpvAccess(ReductionMgr_instance) = this;
    } else {
      DebugM(1, "ReductionMgr::ReductionMgr() - another instance exists!\n");
    }

    // fill in the spanning tree fields
    if (CkMyPe() == 0) {
      myParent = -1;
    } else {
      myParent = (CkMyPe()-1)/REDUCTION_MAX_CHILDREN;
    }
    firstChild = CkMyPe()*REDUCTION_MAX_CHILDREN + 1;
    if (firstChild > CkNumPes()) firstChild = CkNumPes();
    lastChild = firstChild + REDUCTION_MAX_CHILDREN;
    if (lastChild > CkNumPes()) lastChild = CkNumPes();

    // initialize data
    for(int i=0; i<REDUCTION_MAX_SET_ID; i++) {
      reductionSets[i] = 0;
    }

    DebugM(1,"ReductionMgr() instantiated.\n");
}

// destructor
ReductionMgr::~ReductionMgr() {
    for(int i=0; i<REDUCTION_MAX_SET_ID; i++) {
      delete reductionSets[i];
    }

}

// possibly create and return reduction set
ReductionSet* ReductionMgr::getSet(int setID) {
  if ( reductionSets[setID] == 0 ) {
    reductionSets[setID] = new ReductionSet(setID);
    if ( ! isRoot() ) {
      ReductionRegisterMsg *msg = new ReductionRegisterMsg;
      msg->reductionSetID = setID;
      msg->sourceNode = CkMyPe();
      CProxy_ReductionMgr(thisgroup).remoteRegister(msg,myParent);
    }
  }
  return reductionSets[setID];
}

// possibly delete reduction set
void ReductionMgr::delSet(int setID) {
  ReductionSet *set = reductionSets[setID];
  if ( set && ! set->eventsRegistered ) {
    if ( ! isRoot() ) {
      ReductionRegisterMsg *msg = new ReductionRegisterMsg;
      msg->reductionSetID = setID;
      msg->sourceNode = CkMyPe();
      CProxy_ReductionMgr(thisgroup).remoteUnregister(msg,myParent);
    }
    delete set;
    reductionSets[setID] = 0;
  }
}

// register local submit
SubmitReduction* ReductionMgr::willSubmit(int setID) {
  ReductionSet *set = getSet(setID);
  if ( set->dataQueue ) {
    NAMD_die("ReductionMgr::willSubmit called while reductions outstanding!");
  }

  set->eventsRegistered++;

  SubmitReduction *handle = new SubmitReduction;
  handle->reductionSetID = setID;
  handle->sequenceNumber = set->nextSequenceNumber;
  handle->master = this;

  return handle;
}

// unregister local submit
void ReductionMgr::remove(SubmitReduction* handle) {
  int setID = handle->reductionSetID;
  ReductionSet *set = reductionSets[setID];
  if ( set->dataQueue ) {
    NAMD_die("SubmitReduction deleted while reductions outstanding!");
  }

  set->eventsRegistered--;

  delSet(setID);
}

// local submit
void ReductionMgr::submit(SubmitReduction* handle) {
  int setID = handle->reductionSetID;
  ReductionSet *set = reductionSets[setID];
  int seqNum = handle->sequenceNumber;
  int size = handle->dataSize;
  BigReal *data = handle->data;

  mergeAndDeliver(set,seqNum,data,size);
}

// register submit from child
void ReductionMgr::remoteRegister(ReductionRegisterMsg *msg) {

  int setID = msg->reductionSetID;
  ReductionSet *set = getSet(setID);
  if ( set->dataQueue ) {
    NAMD_die("ReductionMgr::willSubmit called while reductions outstanding on parent!");
  }

  set->eventsRegistered++;
  set->addToRemoteSequenceNumber[msg->sourceNode - firstChild]
					= set->nextSequenceNumber;
  delete msg;
}

// unregister submit from child
void ReductionMgr::remoteUnregister(ReductionRegisterMsg *msg) {

  int setID = msg->reductionSetID;
  ReductionSet *set = reductionSets[setID];
  if ( set->dataQueue ) {
    NAMD_die("SubmitReduction deleted while reductions outstanding on parent!");
  }

  set->eventsRegistered--;

  delSet(setID);
  delete msg;
}

// data submitted from child
void ReductionMgr::remoteSubmit(ReductionSubmitMsg *msg) {
  int setID = msg->reductionSetID;
  ReductionSet *set = reductionSets[setID];
  int seqNum = msg->sequenceNumber
	+ set->addToRemoteSequenceNumber[msg->sourceNode - firstChild];
  int size = msg->dataSize;
  BigReal *data = msg->data;

  mergeAndDeliver(set,seqNum,data,size);
  delete msg;
}

// common code for submission and delivery
void ReductionMgr::mergeAndDeliver(
	ReductionSet *set, int seqNum, const BigReal *newData, int size) {
  ReductionSetData *data = set->getData(seqNum);

  // merge in this submission
  data->resize(size);  // extend as needed
  for ( int i = 0; i < size; ++i ) {
    data->data[i] += newData[i];
  }
  data->eventsRemaining--;

  // deliver if all submissions are in
  if ( data->eventsRemaining == set->requireRegistered ) {
    if ( isRoot() ) {
      if ( set->requireRegistered ) {
	if ( set->threadIsWaiting && set->waitingForSequenceNumber == seqNum) {
	  // awaken the thread so it can take the data
	  CthAwaken(set->waitingThread);
	}
      } else {
	NAMD_die("ReductionSet::deliver will never deliver data");
      }
    } else {
      // send data to parent
      int size = data->dataSize;
      ReductionSubmitMsg *msg = new(&size,1) ReductionSubmitMsg;
      msg->reductionSetID = set->reductionSetID;
      msg->sourceNode = CkMyPe();
      msg->sequenceNumber = seqNum;
      msg->dataSize = data->dataSize;
      for ( int i = 0; i < msg->dataSize; ++i ) {
        msg->data[i] = data->data[i];
      }
      CProxy_ReductionMgr(thisgroup).remoteSubmit(msg,myParent);
    }
    set->delData(seqNum);
  }
}

// register require
RequireReduction* ReductionMgr::willRequire(int setID) {
  ReductionSet *set = getSet(setID);
  set->eventsRegistered++;
  set->requireRegistered++;
  if ( set->dataQueue ) {
    NAMD_die("ReductionMgr::willRequire called while reductions outstanding!");
  }

  RequireReduction *handle = new RequireReduction;
  handle->reductionSetID = setID;
  handle->sequenceNumber = set->nextSequenceNumber;
  handle->master = this;

  return handle;
}

// unregister require
void ReductionMgr::remove(RequireReduction* handle) {
  int setID = handle->reductionSetID;
  ReductionSet *set = reductionSets[setID];
  if ( set->dataQueue ) {
    NAMD_die("RequireReduction deleted while reductions outstanding!");
  }

  set->eventsRegistered--;
  set->requireRegistered--;

  delSet(setID);
}

// require the data from a thread
void ReductionMgr::require(RequireReduction* handle) {
  int setID = handle->reductionSetID;
  ReductionSet *set = reductionSets[setID];
  int seqNum = handle->sequenceNumber;
  ReductionSetData *data = set->getData(seqNum);
  if ( data->eventsRemaining > set->requireRegistered ) {
    set->threadIsWaiting = 1;
    set->waitingForSequenceNumber = seqNum;
    set->waitingThread = CthSelf();
    CthSuspend();
  }
  set->threadIsWaiting = 0;
  data->eventsRemaining--;

  if ( handle->dataSize < data->dataSize ) {
    delete [] handle->data;
    handle->data = new BigReal[data->dataSize];
    handle->dataSize = data->dataSize;
  }
  for ( int i = 0; i < data->dataSize; ++i ) {
    handle->data[i] = data->data[i];
  }

  set->delData(seqNum);
}


#include "ReductionMgr.def.h"
// nothing should be placed below here


/***************************************************************************
 * RCS INFORMATION:
 *
 *	$RCSfile $
 *	$Author $	$Locker:  $		$State: Exp $
 *	$Revision: 1.1035 $	$Date: 1999/06/18 17:04:13 $
 *
 ***************************************************************************
 * REVISION HISTORY:
 *
 * $Log: ReductionMgr.C,v $
 * Revision 1.1035  1999/06/18 17:04:13  jim
 * Full dynamics allocation and varsize messages for ReductionMgr.
 *
 * Revision 1.1034  1999/06/17 16:02:54  jim
 * Added error checking for barrier violations.
 *
 * Revision 1.1033  1999/06/17 15:46:16  jim
 * Completely rewrote reduction system to eliminate need for sequence numbers.
 *
 * Revision 1.1032  1999/05/11 23:56:47  brunner
 * Changes for new charm version
 *
 * Revision 1.1031  1999/04/29 15:39:23  jim
 * Added check for extra reduction submissions.
 *
 * Revision 1.1030  1999/02/17 05:43:13  jim
 * Fixed memory leak in more nodes than patches code.
 *
 * Revision 1.1029  1999/02/17 05:29:01  jim
 * Fixed bug in more nodes than patches code.
 *
 * Revision 1.1028  1998/11/30 04:15:27  krishnan
 * Fixed the numNodes > nPatches bug. Dummy reduction is submitted on the nodes with no patches.
 *
 * Revision 1.1027  1998/10/24 19:57:56  jim
 * Eliminated warnings generated by g++ -Wall.
 *
 * Revision 1.1026  1998/05/15 16:19:05  jim
 * Made Controller suspend during load balancing (for reduction system).
 *
 * Revision 1.1025  1998/03/26 23:28:34  jim
 * Small changes for KCC port.  Altered use of strstream in ComputeFreeEnergy.
 *
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
