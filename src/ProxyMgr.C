/**
***  Copyright (c) 1995, 1996, 1997, 1998, 1999, 2000 by
***  The Board of Trustees of the University of Illinois.
***  All rights reserved.
**/

#include "charm++.h"

#include "main.h"
#include "BOCgroup.h"
#include "ProxyMgr.decl.h"
#include "ProxyMgr.h"
#include "Namd.h"
#include "PatchMap.inl"
#include "ProxyPatch.h"
#include "ComputeMap.h"
#include "HomePatch.h"
#include <string.h>
#include "ProcessorPrivate.h"
#include "packmsg.h"

//#define DEBUGM
#define MIN_DEBUG_LEVEL 4
#include "Debug.h"


PACK_MSG(ProxyAtomsMsg,
  PACK(patch);
  PACK_RESIZE(atomIDList);
)
  

PACK_MSG(ProxyDataMsg,
  PACK(patch);
  PACK(flags);
  PACK_RESIZE(positionList);
  if (packmsg_msg->flags.doMolly) PACK_RESIZE(avgPositionList);
)


PACK_MSG(ProxyAllMsg,
  PACK(patch);
  PACK(flags);
  PACK_RESIZE(atomIDList);
  PACK_RESIZE(positionList);
  if (packmsg_msg->flags.doMolly) PACK_RESIZE(avgPositionList);
)


PACK_MSG(ProxyResultMsg,
  PACK(node);
  PACK(patch);
  for ( int j = 0; j < Results::maxNumForces; ++j ) {
    PACK_RESIZE(forceList[j]);
  }
)


ProxyMgr::ProxyMgr() { 
  if (CpvAccess(ProxyMgr_instance)) {
    Namd::die();
  }
  CpvAccess(ProxyMgr_instance) = this;
}

ProxyMgr::~ProxyMgr() { 
  removeProxies();
  CpvAccess(ProxyMgr_instance) = NULL;
}

void ProxyMgr::removeProxies(void)
{
  ProxySetIter pi(proxySet);
  for ( pi = pi.begin(); pi != pi.end(); pi++)
  {
    delete pi->proxyPatch;
  }
  proxySet.clear();
}

// Figure out which proxies we need and create them
void ProxyMgr::createProxies(void)
{
  // Delete the old proxies.
  removeProxies();

  PatchMap *patchMap = PatchMap::Object();
  int numPatches = patchMap->numPatches();
  int myNode = CkMyPe();
  enum PatchFlag { Unknown, Home, NeedProxy };
  int *patchFlag = new int[numPatches]; 
  int i, j;

  // Note all home patches.
  for ( i = 0; i < numPatches; ++i )
  {
    patchFlag[i] = ( patchMap->node(i) == myNode ) ? Home : Unknown;
  }

  // Add all upstream neighbors.
  PatchID neighbors[PatchMap::MaxOneAway];
  for ( i = 0; i < numPatches; ++i )
  {
    if ( patchMap->node(i) != myNode ) 
      continue;
    int numNeighbors = patchMap->upstreamNeighbors(i,neighbors);
    for ( j = 0; j < numNeighbors; ++j )
    {
      if ( ! patchFlag[neighbors[j]] ) {
	patchFlag[neighbors[j]] = NeedProxy;
      }
    }
  }

  // Check all patch-based compute objects.
  ComputeMap *computeMap = ComputeMap::Object();
  int nc = computeMap->numComputes();
  for ( i = 0; i < nc; ++i )
  {
    if ( computeMap->node(i) != myNode || !computeMap->isPatchBased(i) ) 
      continue;
    int numPid = computeMap->numPids(i);
    for ( j = 0; j < numPid; ++j )
    {
      int pid = computeMap->pid(i,j);
      if ( ! patchFlag[pid] ) {
	patchFlag[pid] = NeedProxy;
      }
    }
  }
  
  // Create proxy list
  for ( i = 0; i < numPatches; ++i ) {
    if ( patchFlag[i] == NeedProxy )
    { // create proxy patch
      ProxyPatch *proxy = new ProxyPatch(i);
      proxySet.add(ProxyElem(i, proxy));
      patchMap->registerPatch(i, proxy);
    }
  }
  delete[] patchFlag;
}

void
ProxyMgr::createProxy(PatchID pid) {
  Patch *p = PatchMap::Object()->patch(pid);
  if (!p) {
     DebugM(4,"createProxy("<<pid<<")\n");
     ProxyPatch *proxy = new ProxyPatch(pid);
     proxySet.add(ProxyElem(pid,proxy));
     PatchMap::Object()->registerPatch(pid,proxy);
  }
  else {
     DebugM(4,"createProxy("<<pid<<") found " << p->getPatchID() << "\n");
  }
    
}

void
ProxyMgr::removeProxy(PatchID pid) {
  ProxyElem *p = proxySet.find(ProxyElem(pid));
  if (p) { 
    delete p->proxyPatch;
    proxySet.del(ProxyElem(pid));
  }
}
  
void
ProxyMgr::registerProxy(PatchID pid) {
  // determine which node gets message
  NodeID node = PatchMap::Object()->node(pid);

  RegisterProxyMsg *msg = new RegisterProxyMsg;

  msg->node=CkMyPe();
  msg->patch = pid;

  CProxy_ProxyMgr cp(CpvAccess(BOCclass_group).proxyMgr);
  cp.recvRegisterProxy(msg,node);
}

void
ProxyMgr::recvRegisterProxy(RegisterProxyMsg *msg) {
  HomePatch *homePatch = PatchMap::Object()->homePatch(msg->patch);
  homePatch->registerProxy(msg); // message deleted in registerProxy()
}

void
ProxyMgr::sendResults(ProxyResultMsg *msg) {
  NodeID node = PatchMap::Object()->node(msg->patch);
  CProxy_ProxyMgr cp(CpvAccess(BOCclass_group).proxyMgr);
  cp.recvResults(msg, node);
}

void
ProxyMgr::recvResults(ProxyResultMsg *msg) {
  HomePatch *home = PatchMap::Object()->homePatch(msg->patch);
  home->receiveResults(msg); // delete done in HomePatch::receiveResults()
}

void
ProxyMgr::sendProxyData(ProxyDataMsg *msg, NodeID node) {
  CProxy_ProxyMgr cp(CpvAccess(BOCclass_group).proxyMgr);
  cp.recvProxyData(msg,node);
}

void
ProxyMgr::recvProxyData(ProxyDataMsg *msg) {
  ProxyPatch *proxy = (ProxyPatch *) PatchMap::Object()->patch(msg->patch);
  proxy->receiveData(msg); // deleted in ProxyPatch::receiveAtoms()
}

void
ProxyMgr::sendProxyAtoms(ProxyAtomsMsg *msg, NodeID node) {
  CProxy_ProxyMgr cp(CpvAccess(BOCclass_group).proxyMgr);
  cp.recvProxyAtoms(msg,node);
}

void
ProxyMgr::recvProxyAtoms(ProxyAtomsMsg *msg) {
  ProxyPatch *proxy = (ProxyPatch *) PatchMap::Object()->patch(msg->patch);
  proxy->receiveAtoms(msg); // deleted in ProxyPatch::receiveAtoms()
}

void
ProxyMgr::sendProxyAll(ProxyAllMsg *msg, NodeID node) {
  CProxy_ProxyMgr cp(CpvAccess(BOCclass_group).proxyMgr);
  cp.recvProxyAll(msg,node);
}

void
ProxyMgr::recvProxyAll(ProxyAllMsg *msg) {
  ProxyPatch *proxy = (ProxyPatch *) PatchMap::Object()->patch(msg->patch);
  proxy->receiveAll(msg); // delete done in ProxyPatch::receiveAll()
}

#include "ProxyMgr.def.h"

