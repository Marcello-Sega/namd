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
#include "BroadcastMgr.h"
#include "BroadcastClient.h"
#define MIN_DEBUG_LEVEL 3
// #define DEBUGM
#include "Debug.h"

BroadcastClient::BroadcastClient(int id) {
  this->id = id;
  BroadcastMgr::Object()->subscribe(*this);
  waitForTag = -1;
  suspended = 0;
}

BroadcastClient::~BroadcastClient() {
  BroadcastMgr::Object()->unsubscribe(*this);
}

void 
BroadcastClient::awaken(int theid, int tag) {
  DebugM(1, "awaken() client id = " << id << " tag = " << tag << "\n");
  if (suspended && theid == this->id && tag == waitForTag) {
    CthAwaken(thread);
    suspended = 0; waitForTag = -1;
  }
}

void 
BroadcastClient::suspendFor(int tag) {
  DebugM(1, "suspending() client id = " << id << " tag = " << tag << "\n");
  suspended = 1;
  waitForTag = tag;
  thread = CthSelf();
  CthSuspend();
}

