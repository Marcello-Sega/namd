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

#include "ckdefs.h"
#include "chare.h"
#include "c++interface.h"
#include "BroadcastMgr.h"
#include "BroadcastClient.h"
#define MIN_DEBUG_LEVEL 3
#define DEBUGM
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



/***************************************************************************
 * RCS INFORMATION:
 *
 *	$RCSfile $
 *	$Author $	$Locker:  $		$State: Exp $
 *	$Revision: 1.1 $	$Date: 1997/03/19 11:53:51 $
 *
 ***************************************************************************
 * REVISION HISTORY:
 *
 * $Log: BroadcastClient.C,v $
 * Revision 1.1  1997/03/19 11:53:51  ari
 * Add Broadcast mechanism.
 * Fixed RCS Log entries on files that did not have Log entries.
 * Added some register variables to Molecule and ComputeNonbondedExcl.C
 *
 ***************************************************************************/
