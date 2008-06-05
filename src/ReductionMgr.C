/**
***  Copyright (c) 1995, 1996, 1997, 1998, 1999, 2000 by
***  The Board of Trustees of the University of Illinois.
***  All rights reserved.
**/

/*
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
*/

#include <stdlib.h>
#include <stdio.h>

#include "InfoStream.h"
#include "PatchMap.h"	// for patchMap

#include "Node.h"
#include "SimParameters.h"

#include "ReductionMgr.decl.h"
#include "ReductionMgr.h"

// #define DEBUGM
#define MIN_DEBUG_LEVEL 4
#include "Debug.h"

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
  BigReal data[1];

  static void *alloc(int msgnum, int size, int *array, int priobits) {
    int totalsize = size + array[0]*sizeof(BigReal);
    ReductionSubmitMsg *newMsg = (ReductionSubmitMsg*)
				CkAllocMsg(msgnum,totalsize,priobits);
    // newMsg->data = (BigReal*) ((char*)newMsg + size);
    return (void*)newMsg;
  }

  static void *pack(ReductionSubmitMsg *in) {
    // in->data = (BigReal*) ((char*)in->data - (char*)&(in->data));
    return (void*)in;
  }

  static ReductionSubmitMsg *unpack(void *in) {
    ReductionSubmitMsg *me = (ReductionSubmitMsg*)in;
    // me->data = (BigReal*) ((char*)&(me->data) + (size_t)(me->data));
    return me;
  }

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
#if CHARM_VERSION > 050402
      CProxy_ReductionMgr reductionProxy(thisgroup);
      reductionProxy[myParent].remoteRegister(msg);
#else
      CProxy_ReductionMgr(thisgroup).remoteRegister(msg,myParent);
#endif
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
#if CHARM_VERSION > 050402
      CProxy_ReductionMgr reductionProxy(thisgroup);
      reductionProxy[myParent].remoteUnregister(msg);
#else
      CProxy_ReductionMgr(thisgroup).remoteUnregister(msg,myParent);
#endif
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
  BigReal *curData = data->data;
#ifdef ARCH_POWERPC
#pragma disjoint (*curData,  *newData)
#pragma unroll(4)
#endif
  for ( int i = 0; i < size; ++i ) {
    curData[i] += newData[i];
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
#if CHARM_VERSION > 050402
      CProxy_ReductionMgr reductionProxy(thisgroup);
      reductionProxy[myParent].remoteSubmit(msg);
#else
      CProxy_ReductionMgr(thisgroup).remoteSubmit(msg,myParent);
#endif
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

