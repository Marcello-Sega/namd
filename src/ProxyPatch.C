/**
***  Copyright (c) 1995, 1996, 1997, 1998, 1999, 2000 by
***  The Board of Trustees of the University of Illinois.
***  All rights reserved.
**/

#include "InfoStream.h"

#ifdef USE_COMM_LIB
#include "ComlibManager.h"
#endif

#include "main.decl.h"
#include "main.h"
#include "ProxyPatch.h"
#include "ProxyMgr.decl.h"
#include "ProxyMgr.h"
#include "AtomMap.h"
#include "PatchMap.h"
#include "Priorities.h"

#define MIN_DEBUG_LEVEL 4
//#define  DEBUGM
#include "Debug.h"

ProxyPatch::ProxyPatch(PatchID pd) : 
  Patch(pd), msgBuffer(NULL), msgAllBuffer(NULL)
{
  DebugM(4, "ProxyPatch(" << pd << ") at " << this << "\n");
  ProxyMgr::Object()->registerProxy(patchID);
  numAtoms = -1;
  parent = -1;
  nChild = 0;
  child = new int[proxySpanDim];

#if CMK_PERSISTENT_COMM
  localphs = 0;
  localphs = CmiCreatePersistent(PatchMap::Object()->node(patchID), 300000);
#endif

  // DMK - Atom Separation (water vs. non-water)
  #if NAMD_SeparateWaters != 0
    numWaterAtoms = -1;
  #endif
}

ProxyPatch::~ProxyPatch()
{
  DebugM(4, "ProxyPatch(" << pd << ") deleted at " << this << "\n");
  ProxyMgr::Object()->unregisterProxy(patchID);
  AtomMap::Object()->unregisterIDs(patchID,p.begin(),p.end());
  delete [] child;
#if CMK_PERSISTENT_COMM
  CmiDestoryPersistent(localphs);
  localphs = 0;
#endif
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
  numAtoms = msg->atomIDList.size();
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
  if ( numAtoms == -1 ) { // for new proxies since receiveAtoms is not called
    numAtoms = p.size();

    // DMK - Atom Separation (water vs. non-water)
    #if NAMD_SeparateWaters != 0
      numWaterAtoms = msg->numWaterAtoms;
    #endif

    positionsReady(1);
  } else {
    positionsReady(0);
  }
}

void ProxyPatch::receiveAll(ProxyAllMsg *msg)
{
  DebugM(3, "receiveAll(" << patchID << ")\n");
  if ( boxesOpen )
  {
    // store message in queue (only need one element, though)
    msgAllBuffer = msg;
    return;
  }
  msgAllBuffer = NULL;

  AtomMap::Object()->unregisterIDs(patchID,p.begin(),p.end());
  flags = msg->flags;
  p = msg->positionList;
  numAtoms = p.size();
  p_avg = msg->avgPositionList;

  // DMK - Atom Separation (water vs. non-water)
  #if NAMD_SeparateWaters != 0
    numWaterAtoms = msg->numWaterAtoms;
  #endif

#ifdef MEM_OPT_VERSION
  pExt = msg->extInfoList;
#endif

  delete msg;

  positionsReady(1);
}

void ProxyPatch::sendResults(void)
{
  DebugM(3, "sendResults(" << patchID << ")\n");
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

#if CMK_PERSISTENT_COMM
//  CmiUsePersistentHandle(&localphs, 1);
#endif
  if (proxyRecvSpanning == 0) {
    ProxyResultMsg *msg = new (PRIORITY_SIZE) ProxyResultMsg;
    SET_PRIORITY(msg,flags.sequence,
		PROXY_RESULTS_PRIORITY + PATCH_PRIORITY(patchID));
    msg->node = CkMyPe();
    msg->patch = patchID;
    for ( i = 0; i < Results::maxNumForces; ++i ) 
      msg->forceList[i] = f[i];
    ProxyMgr::Object()->sendResults(msg);
  }
  else {
    ProxyCombinedResultMsg *msg = new (PRIORITY_SIZE) ProxyCombinedResultMsg;
    SET_PRIORITY(msg,flags.sequence,
		PROXY_RESULTS_PRIORITY + PATCH_PRIORITY(patchID));
    msg->nodes.add(CkMyPe());
    msg->patch = patchID;
    for ( i = 0; i < Results::maxNumForces; ++i ) 
      msg->forceList[i] = f[i];
    ProxyMgr::Object()->sendResults(msg);
  }
#if CMK_PERSISTENT_COMM
  CmiUsePersistentHandle(NULL, 0);
#endif
}

void ProxyPatch::setSpanningTree(int p, int *c, int n) { 
  parent=p; nChild = n; nWait = 0;
  for (int i=0; i<n; i++) child[i] = c[i];
//CkPrintf("setSpanningTree: [%d:%d] %d %d:%d %d\n", CkMyPe(), patchID, parent, nChild, child[0], child[1]);
}

int ProxyPatch::getSpanningTreeChild(int *c) { 
  for (int i=0; i<nChild; i++) c[i] = child[i];
  return nChild;
}

ProxyCombinedResultMsg *ProxyPatch::depositCombinedResultMsg(ProxyCombinedResultMsg *msg) {
  nWait++;
  if (nWait == 1) msgCBuffer = msg;
  else {
    NodeIDList::iterator n_i, n_e;
    n_i = msg->nodes.begin();
    n_e = msg->nodes.end();
    for (; n_i!=n_e; ++n_i) msgCBuffer->nodes.add(*n_i);
    for ( int k = 0; k < Results::maxNumForces; ++k )
    {
    register ForceList::iterator r_i;
    r_i = msgCBuffer->forceList[k].begin();
    register ForceList::iterator f_i, f_e;
    f_i = msg->forceList[k].begin();
    f_e = msg->forceList[k].end();
    //    for ( ; f_i != f_e; ++f_i, ++r_i ) *r_i += *f_i;

    int nf = f_e - f_i;
#ifdef ARCH_POWERPC
#pragma disjoint (*f_i, *r_i)
#pragma unroll(4)
#endif
    for (int count = 0; count < nf; count++) {
      r_i[count].x += f_i[count].x;      
      r_i[count].y += f_i[count].y;      
      r_i[count].z += f_i[count].z;
    }

    }
    delete msg;
  }
//CkPrintf("[%d:%d] wait: %d of %d (%d %d %d)\n", CkMyPe(), patchID, nWait, nChild+1, parent, child[0],child[1]);
  if (nWait == nChild + 1) {
    nWait = 0;
    return msgCBuffer;
  }
  return NULL;
}

