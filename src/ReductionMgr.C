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
  int dataSize;
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

ReductionSet::ReductionSet(int setID, int size) {
  if ( setID == REDUCTIONS_BASIC ) {
    if ( size != -1 ) {
      NAMD_bug("ReductionSet size specified for REDUCTIONS_BASIC.");
    }
    size = REDUCTION_MAX_RESERVED;
  }
  if ( size == -1 ) NAMD_bug("ReductionSet size not specified.");
  dataSize = size;
  reductionSetID = setID;
  nextSequenceNumber = 0;
  submitsRegistered = 0;
  dataQueue = 0;
  requireRegistered = 0;
  threadIsWaiting = 0;
}

ReductionSet::~ReductionSet() {

  ReductionSetData *current = dataQueue;

  while ( current ) {
    ReductionSetData *next = current->next;
    delete current;
    current = next;
  }
}

// possibly create and return data for a particular seqNum
ReductionSetData* ReductionSet::getData(int seqNum) {

  ReductionSetData **current = &dataQueue;

  while ( *current ) {
    if ( (*current)->sequenceNumber == seqNum ) return *current;
    current = &((*current)->next);
  }

//iout << "seq " << seqNum << " created on " << CkMyPe() << "\n" << endi;
  *current = new ReductionSetData(seqNum, dataSize);
  return *current;
}

// possibly delete data for a particular seqNum
ReductionSetData* ReductionSet::removeData(int seqNum) {

  ReductionSetData **current = &dataQueue;

  while ( *current ) {
    if ( (*current)->sequenceNumber == seqNum ) break;
    current = &((*current)->next);
  }

  if ( ! *current ) { NAMD_die("ReductionSet::removeData on missing seqNum"); }

  ReductionSetData *toremove = *current;
  *current = (*current)->next;
  return toremove;
}

// constructor
ReductionMgr::ReductionMgr() {
    if (CkpvAccess(ReductionMgr_instance) == 0) {
      CkpvAccess(ReductionMgr_instance) = this;
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
ReductionSet* ReductionMgr::getSet(int setID, int size) {
  if ( reductionSets[setID] == 0 ) {
    reductionSets[setID] = new ReductionSet(setID,size);
    if ( ! isRoot() ) {
      ReductionRegisterMsg *msg = new ReductionRegisterMsg;
      msg->reductionSetID = setID;
      msg->dataSize = size;
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
  if ( set && ! set->submitsRegistered & ! set->requireRegistered ) {
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
SubmitReduction* ReductionMgr::willSubmit(int setID, int size) {
  ReductionSet *set = getSet(setID, size);
  ReductionSetData *data = set->getData(set->nextSequenceNumber);
  if ( data->submitsRecorded ) {
    NAMD_die("ReductionMgr::willSubmit called while reductions outstanding!");
  }

  set->submitsRegistered++;

  SubmitReduction *handle = new SubmitReduction;
  handle->reductionSetID = setID;
  handle->sequenceNumber = set->nextSequenceNumber;
  handle->master = this;
  handle->data = data->data;

  return handle;
}

// unregister local submit
void ReductionMgr::remove(SubmitReduction* handle) {
  int setID = handle->reductionSetID;
  ReductionSet *set = reductionSets[setID];
  if ( set->getData(set->nextSequenceNumber)->submitsRecorded ) {
    NAMD_die("SubmitReduction deleted while reductions outstanding!");
  }

  set->submitsRegistered--;

  delSet(setID);
}

// local submit
void ReductionMgr::submit(SubmitReduction* handle) {
  int setID = handle->reductionSetID;
  int seqNum = handle->sequenceNumber;
  ReductionSet *set = reductionSets[setID];
  ReductionSetData *data = set->getData(seqNum);

  data->submitsRecorded++;
  if ( data->submitsRecorded == set->submitsRegistered ) {
    mergeAndDeliver(set,seqNum);
  }

  handle->sequenceNumber = ++seqNum;
  handle->data = set->getData(seqNum)->data;
}

// register submit from child
void ReductionMgr::remoteRegister(ReductionRegisterMsg *msg) {

  int setID = msg->reductionSetID;
  int size = msg->dataSize;
  ReductionSet *set = getSet(setID,size);
  if ( set->getData(set->nextSequenceNumber)->submitsRecorded ) {
    NAMD_die("ReductionMgr::remoteRegister called while reductions outstanding on parent!");
  }

  set->submitsRegistered++;
  set->addToRemoteSequenceNumber[msg->sourceNode - firstChild]
					= set->nextSequenceNumber;
  delete msg;
}

// unregister submit from child
void ReductionMgr::remoteUnregister(ReductionRegisterMsg *msg) {

  int setID = msg->reductionSetID;
  ReductionSet *set = reductionSets[setID];
  if ( set->getData(set->nextSequenceNumber)->submitsRecorded ) {
    NAMD_die("SubmitReduction deleted while reductions outstanding on parent!");
  }

  set->submitsRegistered--;

  delSet(setID);
  delete msg;
}

// data submitted from child
void ReductionMgr::remoteSubmit(ReductionSubmitMsg *msg) {
  int setID = msg->reductionSetID;
  ReductionSet *set = reductionSets[setID];
  int seqNum = msg->sequenceNumber
	+ set->addToRemoteSequenceNumber[msg->sourceNode - firstChild];

//iout << "seq " << seqNum << " from " << msg->sourceNode << " received on " << CkMyPe() << "\n" << endi;
  int size = msg->dataSize;
  if ( size != set->dataSize ) {
    NAMD_bug("ReductionMgr::remoteSubmit data sizes do not match.");
  }

  BigReal *newData = msg->data;
  ReductionSetData *data = set->getData(seqNum);
  BigReal *curData = data->data;
#ifdef ARCH_POWERPC
#pragma disjoint (*curData,  *newData)
#pragma unroll(4)
#endif
  for ( int i = 0; i < size; ++i ) {
    curData[i] += newData[i];
  }
  delete msg;

  data->submitsRecorded++;
  if ( data->submitsRecorded == set->submitsRegistered ) {
    mergeAndDeliver(set,seqNum);
  }
}

// common code for submission and delivery
void ReductionMgr::mergeAndDeliver(ReductionSet *set, int seqNum) {

//iout << "seq " << seqNum << " complete on " << CkMyPe() << "\n" << endi;
 
    set->nextSequenceNumber++; // should match all clients

    ReductionSetData *data = set->getData(seqNum);
    if ( data->submitsRecorded != set->submitsRegistered ) {
      NAMD_bug("ReductionMgr::mergeAndDeliver not ready to deliver.");
    }

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
      int size = set->dataSize;
      ReductionSubmitMsg *msg = new(&size,1) ReductionSubmitMsg;
      msg->reductionSetID = set->reductionSetID;
      msg->sourceNode = CkMyPe();
      msg->sequenceNumber = seqNum;
      msg->dataSize = set->dataSize;
      for ( int i = 0; i < msg->dataSize; ++i ) {
        msg->data[i] = data->data[i];
      }
#if CHARM_VERSION > 050402
      CProxy_ReductionMgr reductionProxy(thisgroup);
      reductionProxy[myParent].remoteSubmit(msg);
#else
      CProxy_ReductionMgr(thisgroup).remoteSubmit(msg,myParent);
#endif
      delete set->removeData(seqNum);
    }

}

// register require
RequireReduction* ReductionMgr::willRequire(int setID, int size) {
  ReductionSet *set = getSet(setID,size);
  set->requireRegistered++;
  if ( set->getData(set->nextSequenceNumber)->submitsRecorded ) {
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
  if ( set->getData(set->nextSequenceNumber)->submitsRecorded ) {
    NAMD_die("RequireReduction deleted while reductions outstanding!");
  }

  set->requireRegistered--;

  delSet(setID);
}

// require the data from a thread
void ReductionMgr::require(RequireReduction* handle) {
  int setID = handle->reductionSetID;
  ReductionSet *set = reductionSets[setID];
  int seqNum = handle->sequenceNumber;
  ReductionSetData *data = set->getData(seqNum);
  if ( data->submitsRecorded < set->submitsRegistered ) {
    set->threadIsWaiting = 1;
    set->waitingForSequenceNumber = seqNum;
    set->waitingThread = CthSelf();
//iout << "seq " << seqNum << " waiting\n" << endi;
    CthSuspend();
  }
  set->threadIsWaiting = 0;

//iout << "seq " << seqNum << " consumed\n" << endi;
  delete handle->currentData;
  handle->currentData = set->removeData(seqNum);
  handle->data = handle->currentData->data;
  handle->sequenceNumber = ++seqNum;
}


#include "ReductionMgr.def.h"
// nothing should be placed below here

