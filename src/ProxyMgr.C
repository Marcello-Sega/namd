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

int proxySendSpanning = 0;
int proxyRecvSpanning = 0;

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
  PACK_RESIZE(positionList);
  if (packmsg_msg->flags.doMolly) PACK_RESIZE(avgPositionList);
)

PACK_MSG(ProxySpanningTreeMsg,
  PACK(patch);
  PACK(node);
  PACK_RESIZE(tree);
)

void* ProxyResultMsg::pack(ProxyResultMsg *msg) {
  int msg_size = 0;
  msg_size += sizeof(msg->node);
  msg_size += sizeof(msg->patch);
  int j;
  for ( j = 0; j < Results::maxNumForces; ++j ) {
    int array_size = msg->forceList[j].size();
    msg_size += sizeof(array_size);
    msg_size += array_size * sizeof(char);
    Force* f = msg->forceList[j].begin();
    int nonzero_count = 0;
    for ( int i = 0; i < array_size; ++i ) {
      if ( f[i].x != 0. || f[i].y != 0. || f[i].z != 0. ) { ++nonzero_count; }
    }
    msg_size += nonzero_count * sizeof(Force);
  }

  void *msg_buf = CkAllocBuffer(msg,msg_size);
  char *msg_cur = (char *)msg_buf;

  memcpy((void*)msg_cur,(void*)(&(msg->node)),sizeof(msg->node));
  msg_cur += sizeof(msg->node);
  memcpy((void*)msg_cur,(void*)(&(msg->patch)),sizeof(msg->patch));
  msg_cur += sizeof(msg->patch);
  for ( j = 0; j < Results::maxNumForces; ++j ) {
    int array_size = msg->forceList[j].size();
    memcpy((void*)msg_cur,(void*)(&array_size),sizeof(array_size));
    msg_cur += sizeof(array_size);
    char *nonzero = msg_cur;
    msg_cur += array_size * sizeof(char);
    Force* f = msg->forceList[j].begin();
    for ( int i = 0; i < array_size; ++i ) {
      if ( f[i].x != 0. || f[i].y != 0. || f[i].z != 0. ) {
        nonzero[i] = 1;
        memcpy((void*)msg_cur,(void*)(f+i),sizeof(Force));
        msg_cur += sizeof(Force);
      } else {
        nonzero[i] = 0;
      }
    }
  }

  delete msg;
  return msg_buf;
}

ProxyResultMsg* ProxyResultMsg::unpack(void *ptr) {
  void *vmsg = CkAllocBuffer(ptr,sizeof(ProxyResultMsg));
  ProxyResultMsg *msg = new (vmsg) ProxyResultMsg;
  char *msg_cur = (char*)ptr;

  memcpy((void*)(&(msg->node)),(void*)msg_cur,sizeof(msg->node));
  msg_cur += sizeof(msg->node);
  memcpy((void*)(&(msg->patch)),(void*)msg_cur,sizeof(msg->patch));
  msg_cur += sizeof(msg->patch);
  int j;
  for ( j = 0; j < Results::maxNumForces; ++j ) {
    int array_size;
    memcpy((void*)(&array_size),(void*)msg_cur,sizeof(array_size));
    msg_cur += sizeof(array_size);
    msg->forceList[j].resize(array_size);
    char *nonzero = msg_cur;
    msg_cur += array_size * sizeof(char);
    Force* f = msg->forceList[j].begin();
    for ( int i = 0; i < array_size; ++i ) {
      if ( nonzero[i] ) {
        memcpy((void*)(f+i),(void*)msg_cur,sizeof(Force));
        msg_cur += sizeof(Force);
      } else {
        f[i].x = 0.;  f[i].y = 0.;  f[i].z = 0.;
      }
    }
  }

  CkFreeMsg(ptr);
  return msg;
}


// for spanning tree
void* ProxyCombinedResultMsg::pack(ProxyCombinedResultMsg *msg) {
  int msg_size = 0;
  msg_size += sizeof(int) + msg->nodes.size()*sizeof(NodeID);
  msg_size += sizeof(msg->patch);
  int j;
  for ( j = 0; j < Results::maxNumForces; ++j ) {
    int array_size = msg->forceList[j].size();
    msg_size += sizeof(array_size);
    msg_size += array_size * sizeof(char);
    Force* f = msg->forceList[j].begin();
    int nonzero_count = 0;
    for ( int i = 0; i < array_size; ++i ) {
      if ( f[i].x != 0. || f[i].y != 0. || f[i].z != 0. ) { ++nonzero_count; }
    }
    msg_size += nonzero_count * sizeof(Force);
  }

  void *msg_buf = CkAllocBuffer(msg,msg_size);
  char *msg_cur = (char *)msg_buf;

  int nodeSize = msg->nodes.size();
  memcpy((void*)msg_cur,(void*)(&nodeSize), sizeof(nodeSize));
  msg_cur += sizeof(nodeSize);
  for (int i=0; i<nodeSize; i++) {
    memcpy((void*)msg_cur,(void*)(&msg->nodes[i]), sizeof(NodeID));
    msg_cur += sizeof(NodeID);
  }
  memcpy((void*)msg_cur,(void*)(&(msg->patch)),sizeof(msg->patch));
  msg_cur += sizeof(msg->patch);
  for ( j = 0; j < Results::maxNumForces; ++j ) {
    int array_size = msg->forceList[j].size();
    memcpy((void*)msg_cur,(void*)(&array_size),sizeof(array_size));
    msg_cur += sizeof(array_size);
    char *nonzero = msg_cur;
    msg_cur += array_size * sizeof(char);
    Force* f = msg->forceList[j].begin();
    for ( int i = 0; i < array_size; ++i ) {
      if ( f[i].x != 0. || f[i].y != 0. || f[i].z != 0. ) {
        nonzero[i] = 1;
        memcpy((void*)msg_cur,(void*)(f+i),sizeof(Force));
        msg_cur += sizeof(Force);
      } else {
        nonzero[i] = 0;
      }
    }
  }

  delete msg;
  return msg_buf;
}

ProxyCombinedResultMsg* ProxyCombinedResultMsg::unpack(void *ptr) {
  void *vmsg = CkAllocBuffer(ptr,sizeof(ProxyCombinedResultMsg));
  ProxyCombinedResultMsg *msg = new (vmsg) ProxyCombinedResultMsg;
  char *msg_cur = (char*)ptr;

  int nodeSize;
  memcpy((void*)(&nodeSize),(void*)msg_cur,sizeof(nodeSize));
  msg_cur += sizeof(nodeSize);
  for (int i=0; i<nodeSize; i++) {
    msg->nodes.add(*(int *)msg_cur);
    msg_cur += sizeof(NodeID);
  }
  memcpy((void*)(&(msg->patch)),(void*)msg_cur,sizeof(msg->patch));
  msg_cur += sizeof(msg->patch);
  int j;
  for ( j = 0; j < Results::maxNumForces; ++j ) {
    int array_size;
    memcpy((void*)(&array_size),(void*)msg_cur,sizeof(array_size));
    msg_cur += sizeof(array_size);
    msg->forceList[j].resize(array_size);
    char *nonzero = msg_cur;
    msg_cur += array_size * sizeof(char);
    Force* f = msg->forceList[j].begin();
    for ( int i = 0; i < array_size; ++i ) {
      if ( nonzero[i] ) {
        memcpy((void*)(f+i),(void*)msg_cur,sizeof(Force));
        msg_cur += sizeof(Force);
      } else {
        f[i].x = 0.;  f[i].y = 0.;  f[i].z = 0.;
      }
    }
  }

  CkFreeMsg(ptr);
  return msg;
}

ProxyMgr::ProxyMgr() { 
  if (CpvAccess(ProxyMgr_instance)) {
    NAMD_bug("Tried to create ProxyMgr twice.");
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

void ProxyMgr::removeUnusedProxies(void)
{
  ResizeArray<PatchID> toDelete;
  ProxySetIter pi(proxySet);
  for ( pi = pi.begin(); pi != pi.end(); pi++)
  {
    if ( pi->proxyPatch->getNumComputes() == 0 ) {
      toDelete.add(pi->patchID);
    }
  }
  PatchID *pidi = toDelete.begin();
  for ( ; pidi != toDelete.end(); ++pidi ) {
    removeProxy(*pidi);
  }
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
    PatchMap::Object()->unregisterPatch(pid,p->proxyPatch);
    delete p->proxyPatch;
    proxySet.del(ProxyElem(pid));
    iout << iINFO << "Removing unused proxy " << pid << " on " << iPE << ".\n" << endi;
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
ProxyMgr::unregisterProxy(PatchID pid) {
  // determine which node gets message
  NodeID node = PatchMap::Object()->node(pid);

  UnregisterProxyMsg *msg = new UnregisterProxyMsg;

  msg->node=CkMyPe();
  msg->patch = pid;

  CProxy_ProxyMgr cp(CpvAccess(BOCclass_group).proxyMgr);
  cp.recvUnregisterProxy(msg,node);
}

void
ProxyMgr::recvUnregisterProxy(UnregisterProxyMsg *msg) {
  HomePatch *homePatch = PatchMap::Object()->homePatch(msg->patch);
  homePatch->unregisterProxy(msg); // message deleted in registerProxy()
}

void 
ProxyMgr::buildProxySpanningTree()
{
  PatchIDList pids;
  PatchMap::Object()->homePatchIDList(pids);
  for (int i=0; i<pids.size(); i++) {
    HomePatch *home = PatchMap::Object()->homePatch(pids[i]);
    if (home == NULL) CkPrintf("ERROR: homepatch NULL\n");
    home->buildSpanningTree();
  }
}

void 
ProxyMgr::sendSpanningTree(ProxySpanningTreeMsg *msg) {
  CProxy_ProxyMgr cp(CpvAccess(BOCclass_group).proxyMgr);
  cp.recvSpanningTree(msg, msg->tree[0]);
}

void 
ProxyMgr::recvSpanningTree(ProxySpanningTreeMsg *msg) {
  int size = msg->tree.size();
  int child[PROXY_SPAN_DIM];
  int nChild = 0;
  int i;
  ProxyPatch *proxy = (ProxyPatch *) PatchMap::Object()->patch(msg->patch);
  for (i=0; i<PROXY_SPAN_DIM; i++) {
    if (size > i+1) { child[i] = msg->tree[i+1]; nChild++; }
  }
  if (!PatchMap::Object()->homePatch(msg->patch)) {
    proxy->setSpanningTree(msg->node, child, nChild);
  }

  // build subtree and pass down
  if (nChild == 0) return;
//CkPrintf("[%d] %d:(%d) %d %d %d %d %d\n", CkMyPe(), msg->patch, size, msg->tree[0], msg->tree[1], msg->tree[2], msg->tree[3], msg->tree[4]);
  CProxy_ProxyMgr cp(CpvAccess(BOCclass_group).proxyMgr);
  NodeIDList tree[PROXY_SPAN_DIM];
  int level = 1, index=1;
  int done = 0;
  while (!done) {
    for (int n=0; n<nChild; n++) {
      if (done) break;
      for (int j=0; j<level; j++) {
       if (index >= size) { done = 1; break; }
       tree[n].add(msg->tree[index]);
       index++;
      }
    }
    level *=PROXY_SPAN_DIM;
  }

  ProxyMgr *proxyMgr = ProxyMgr::Object();
  for (i=0; i<PROXY_SPAN_DIM; i++) {
    if (tree[i].size()) {
      ProxySpanningTreeMsg *cmsg = new ProxySpanningTreeMsg;
      cmsg->patch = msg->patch;
      cmsg->node = CkMyPe();
      cmsg->tree = tree[i];
      proxyMgr->sendSpanningTree(cmsg);
    }
  }

  delete msg;
}

void
ProxyMgr::sendResults(ProxyResultMsg *msg) {
  CProxy_ProxyMgr cp(CpvAccess(BOCclass_group).proxyMgr);
  NodeID node = PatchMap::Object()->node(msg->patch);
  cp.recvResults(msg, node);
}

void
ProxyMgr::recvResults(ProxyResultMsg *msg) {
  HomePatch *home = PatchMap::Object()->homePatch(msg->patch);
  home->receiveResults(msg); // delete done in HomePatch::receiveResults()
}

void
ProxyMgr::sendResults(ProxyCombinedResultMsg *msg) {
  ProxyPatch *patch = (ProxyPatch *)PatchMap::Object()->patch(msg->patch);
  ProxyCombinedResultMsg *cMsg = patch->depositCombinedResultMsg(msg);
  if (cMsg) {
    CProxy_ProxyMgr cp(CpvAccess(BOCclass_group).proxyMgr);
    cp.recvResults(cMsg, patch->getSpanningTreeParent());
  }
}

void
ProxyMgr::recvResults(ProxyCombinedResultMsg *msg) {
  HomePatch *home = PatchMap::Object()->homePatch(msg->patch);
  if (home) {
    home->receiveResults(msg); // delete done in HomePatch::receiveResults()
  }
  else {
    ProxyPatch *patch = (ProxyPatch *)PatchMap::Object()->patch(msg->patch);
    ProxyCombinedResultMsg *cMsg = patch->depositCombinedResultMsg(msg);
    if (cMsg) {
      CProxy_ProxyMgr cp(CpvAccess(BOCclass_group).proxyMgr);
      cp.recvResults(cMsg, patch->getSpanningTreeParent());
    }
  }
}

void
ProxyMgr::sendProxyData(ProxyDataMsg *msg, int pcnt, int *pids) {
  CProxy_ProxyMgr cp(CpvAccess(BOCclass_group).proxyMgr);
  cp.recvProxyData(msg,pcnt,pids);
}

void
ProxyMgr::recvProxyData(ProxyDataMsg *msg) {
  ProxyPatch *proxy = (ProxyPatch *) PatchMap::Object()->patch(msg->patch);
  if (proxySendSpanning == 1) {
    // copy the message and send to spanning children
    int npid = 0;
    int pids[PROXY_SPAN_DIM];
    if (npid = proxy->getSpanningTreeChild(pids)) {
      ProxyDataMsg *newmsg = new(sizeof(int)*8) ProxyDataMsg;
      CkSetQueueing(newmsg, CK_QUEUEING_IFIFO);
      *((int*) CkPriorityPtr(newmsg)) = *((int*) CkPriorityPtr(msg));
      newmsg->patch = msg->patch;
      newmsg->flags = msg->flags;
      newmsg->positionList = msg->positionList;
      newmsg->avgPositionList = msg->avgPositionList;
      ProxyMgr::Object()->sendProxyData(newmsg,npid,pids);
    }
  }
  proxy->receiveData(msg); // deleted in ProxyPatch::receiveAtoms()
}

void
ProxyMgr::sendProxyAll(ProxyAllMsg *msg, int pcnt, int *pids) {
  CProxy_ProxyMgr cp(CpvAccess(BOCclass_group).proxyMgr);
  cp.recvProxyAll(msg,pcnt,pids);
}

void
ProxyMgr::recvProxyAll(ProxyAllMsg *msg) {
  ProxyPatch *proxy = (ProxyPatch *) PatchMap::Object()->patch(msg->patch);
  if (proxySendSpanning == 1) {
    // copy the message and send to spanning children
    int npid = 0;
    int pids[PROXY_SPAN_DIM];
    if (npid = proxy->getSpanningTreeChild(pids)) {
      ProxyAllMsg *newmsg = new(sizeof(int)*8) ProxyAllMsg;
      CkSetQueueing(newmsg, CK_QUEUEING_IFIFO);
      *((int*) CkPriorityPtr(newmsg)) = *((int*) CkPriorityPtr(msg));
      newmsg->patch = msg->patch;
      newmsg->flags = msg->flags;
      newmsg->positionList = msg->positionList;
      newmsg->avgPositionList = msg->avgPositionList;
      ProxyMgr::Object()->sendProxyAll(newmsg,npid,pids);
    }
  }
  proxy->receiveAll(msg); // delete done in ProxyPatch::receiveAll()
}

#include "ProxyMgr.def.h"

