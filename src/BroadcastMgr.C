//-*-c++-*-
/***************************************************************************/
/*              (C) Copyright 1996,1997 The Board of Trustees of the       */
/*                          University of Illinois                         */
/*                           All Rights Reserved                           */
/***************************************************************************/
/***************************************************************************
 * DESCRIPTION: Coordinates broadcast of a data type from a Controller/Seq
 *		to all other Controller/Sequencer type objects (they must
 *		run in a thread!)
 ***************************************************************************/

#include "charm++.h"

#include "UniqueSet.h"
#include "UniqueSetIter.h"
#include "BroadcastMgr.decl.h"
#include "BroadcastMgr.h"
#include "BroadcastClient.h"
#include "BroadcastObject.h"
#include "ProcessorPrivate.h"
#define MIN_DEBUG_LEVEL 3
//#define DEBUGM
#include "Debug.h"

BroadcastMgr::~BroadcastMgr(void) {
  UniqueSetIter<BOID> boidIter(boid);
  for (boidIter = boidIter.begin(); boidIter != boidIter.end(); boidIter++) {
    delete boidIter->broadcastSet;
    if (boidIter->taggedMsg) {
      UniqueSetIter<TaggedMsg> tmIter(*(boidIter->taggedMsg));
      for (tmIter = tmIter.begin(); tmIter != tmIter.end(); tmIter++) {
	delete tmIter->msg;
      }
      delete boidIter->taggedMsg;
    }
  }
}


void *
BroadcastMgr::getbuf(BroadcastClient &b, int tag) {
  void *msg = 0;
  TaggedMsg *tm;
  BOID* boidTmp = boid.find(BOID(b.id));
  if (!boidTmp) {
    return(NULL);
  }
  if ( (tm = (boidTmp->taggedMsg->find(TaggedMsg(tag)))) ) {
    msg = (void *)new char[tm->msgSize];
    memcpy(msg, tm->msg, tm->msgSize);
    if (!--(tm->counter)) {
      (boid.find(BOID(b.id)))->taggedMsg->del(TaggedMsg(tag));
    }
  }
  return(msg);
}


void 
BroadcastMgr::send(BroadcastClient &b, int tag, void *buf, size_t size) {
  BroadcastMsg* msg = new BroadcastMsg;
  msg->msg = buf;
  msg->size = (int)size;
  msg->tag = tag;
  msg->id = b.id;
  msg->node = CkMyPe();
  CProxy_BroadcastMgr(thisgroup).recvBroadcast(msg);
}

void 
BroadcastMgr::subscribe(BroadcastClient &bc) {
  BOID *b;
  if (!(b = boid.find(BOID(bc.id)))) {
    boid.add(BOID(bc.id));
    b = boid.find(BOID(bc.id));
    b->broadcastSet = new UniqueSet<BroadcastClientElem>;
    b->taggedMsg = new UniqueSet<TaggedMsg>;
  }
  b->broadcastSet->add(BroadcastClientElem(&bc));
}

void 
BroadcastMgr::unsubscribe(BroadcastClient &bc) {
  BOID *b;
  if ( (b = boid.find(BOID(bc.id))) ) {
    b->broadcastSet->del(BroadcastClientElem(&bc));
    if (!b->broadcastSet->size()) {
      delete b->broadcastSet;
      b->broadcastSet = 0;
      UniqueSetIter<TaggedMsg> tmIter(*(b->taggedMsg));
      for (tmIter = tmIter.begin(); tmIter != tmIter.end(); tmIter++) {
	delete[] tmIter->msg;
      }
      delete b->taggedMsg;
      b->taggedMsg = 0;
    }
  }
}

void 
BroadcastMgr::recvBroadcast(BroadcastMsg *msg) {
  BOID *b;
  int counter;
  // Check if msg->id has any registrants
  if ( (b = boid.find(BOID(msg->id))) ) {
    // add message to taggedMsg container
    counter = b->broadcastSet->size();
    if (msg->node == CkMyPe()) counter--; // get rid of sender
    b->taggedMsg->add(TaggedMsg(msg->tag,msg->size,counter,msg->msg));

    // inform all registrants of mew message
    UniqueSetIter<BroadcastClientElem> bcIter(*(b->broadcastSet));
    for (bcIter = bcIter.begin(); bcIter != bcIter.end(); bcIter++) {
      bcIter->broadcastClient->awaken(msg->id, msg->tag);
    }
  }
  delete msg;
}

#include "BroadcastMgr.def.h"

