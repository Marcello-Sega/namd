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

#include "charm++.h"

#include "main.decl.h"
#include "main.h"
#include "ProxyPatch.h"
#include "ProxyMgr.decl.h"
#include "ProxyMgr.h"
#include "AtomMap.h"

#define MIN_DEBUG_LEVEL 4
//#define  DEBUGM
#include "Debug.h"

ProxyPatch::ProxyPatch(PatchID pd) : 
  Patch(pd), msgBuffer(NULL), msgAllBuffer(NULL)
{
  DebugM(4, "ProxyPatch(" << pd << ") at " << this << "\n");
  ProxyMgr::Object()->registerProxy(patchID);
}

void ProxyPatch::boxClosed(int box)
{
  if ( box == 1 ) {
    sendResults();
  }
  if ( ! --boxesOpen ) {
    DebugM(2,patchID << ": " << "Checking message buffer.\n");
    if ( msgBuffer ) {
      DebugM(3,"Patch " << patchID << " processing buffered proxy data.\n");
      receiveData(msgBuffer);
    } else if (msgAllBuffer ) {
      DebugM(3,"Patch " << patchID << " processing buffered proxy ALL data.\n");
      receiveAll(msgAllBuffer);
    }
  }
  else {
    DebugM(3,"ProxyPatch " << patchID << ": " << boxesOpen << " boxes left to close.\n");
  }
}

void ProxyPatch::receiveAtoms(ProxyAtomsMsg *msg)
{
  DebugM(3, "receiveAtoms(" << patchID << ")\n");
  loadAtoms(msg->atomIDList);
  AtomMap::Object()->registerIDs(patchID,msg->atomIDList);
  delete msg;
}

void ProxyPatch::receiveData(ProxyDataMsg *msg)
{
  DebugM(3, "receiveData(" << patchID << ")\n");
  if ( boxesOpen )
  {
    // store message in queue (only need one element, though)
    msgBuffer = msg;
    return;
  }
  msgBuffer = NULL;
  flags = msg->flags;
  p = msg->positionList;
  p_avg = msg->avgPositionList;
  delete msg;
  positionsReady(0);
}

void ProxyPatch::receiveAll(ProxyAllMsg *msg)
{
  DebugM(3, "receiveData(" << patchID << ")\n");
  if ( boxesOpen )
  {
    // store message in queue (only need one element, though)
    msgAllBuffer = msg;
    return;
  }
  msgAllBuffer = NULL;

  AtomMap::Object()->unregisterIDs(patchID,atomIDList);
  loadAtoms(msg->atomIDList);
  AtomMap::Object()->registerIDs(patchID,msg->atomIDList);
  flags = msg->flags;
  p = msg->positionList;
  p_avg = msg->avgPositionList;

  delete msg;

  positionsReady(1);
}

void ProxyPatch::sendResults(void)
{
  DebugM(3, "sendResults(" << patchID << ")\n");
  ProxyResultMsg *msg = new ProxyResultMsg;
  msg->node = CkMyPe();
  msg->patch = patchID;
  register int i = 0;
  register ForceList::iterator f_i, f_e, f2_i;
  for ( i = Results::normal + 1 ; i <= flags.maxForceMerged; ++i ) {
    f_i = f[Results::normal].begin(); f_e = f[Results::normal].end();
    f2_i = f[i].begin();
    for ( ; f_i != f_e; ++f_i, ++f2_i ) *f_i += *f2_i;
    f[i].resize(0);
  }
  for ( i = flags.maxForceUsed + 1; i < Results::maxNumForces; ++i )
    f[i].resize(0);
  for ( i = 0; i < Results::maxNumForces; ++i ) 
    msg->forceList[i] = f[i];
  ProxyMgr::Object()->sendResults(msg);
}

