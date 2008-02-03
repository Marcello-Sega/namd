/**
***  Copyright (c) 1995, 1996, 1997, 1998, 1999, 2000 by
***  The Board of Trustees of the University of Illinois.
***  All rights reserved.
**/

#include "InfoStream.h"
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
#include "Priorities.h"

//#define DEBUGM
#define MIN_DEBUG_LEVEL 4
#include "Debug.h"

int proxySendSpanning	= 0;
int proxyRecvSpanning	= 0;
int proxySpanDim	= 9;

PACK_MSG(ProxyAtomsMsg,
  PACK(patch);
  PACK_RESIZE(atomIDList);
)


// DMK - Atom Separation (water vs. non-water)
#if NAMD_SeparateWaters != 0
  
PACK_MSG(ProxyDataMsg,
  PACK(patch);
  PACK(flags);
  PACK_RESIZE(positionList);
  if (packmsg_msg->flags.doMolly) PACK_RESIZE(avgPositionList);
  PACK(numWaterAtoms);
)

#else // NAMD_SeparateWaters == 0

PACK_MSG(ProxyDataMsg,
  PACK(patch);
  PACK(flags);
  PACK_RESIZE(positionList);
  if (packmsg_msg->flags.doMolly) PACK_RESIZE(avgPositionList);
)

#endif


// DMK - Atom Separation (water vs. non-water)
#if NAMD_SeparateWaters != 0

#ifdef MEM_OPT_VERSION
PACK_MSG(ProxyAllMsg,
  PACK(patch);
  PACK(flags);
  PACK_RESIZE(positionList);
  if (packmsg_msg->flags.doMolly) PACK_RESIZE(avgPositionList);
  PACK_RESIZE(extInfoList);
  PACK(numWaterAtoms);
)
#else
PACK_MSG(ProxyAllMsg,
  PACK(patch);
  PACK(flags);
  PACK_RESIZE(positionList);
  if (packmsg_msg->flags.doMolly) PACK_RESIZE(avgPositionList);
  PACK(numWaterAtoms);
)
#endif

#else // NAMD_SeparateWaters == 0

#ifdef MEM_OPT_VERSION
PACK_MSG(ProxyAllMsg,
  PACK(patch);
  PACK(flags);
  PACK_RESIZE(positionList);
  if (packmsg_msg->flags.doMolly) PACK_RESIZE(avgPositionList);
  PACK_RESIZE(extInfoList);
)
#else
PACK_MSG(ProxyAllMsg,
  PACK(patch);
  PACK(flags);
  PACK_RESIZE(positionList);
  if (packmsg_msg->flags.doMolly) PACK_RESIZE(avgPositionList);
)
#endif

#endif


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
    msg_size = ALIGN_8 (msg_size);
    Force* f = msg->forceList[j].begin();
    int nonzero_count = 0;
    for ( int i = 0; i < array_size; ++i ) {
      if ( f[i].x != 0. || f[i].y != 0. || f[i].z != 0. ) { ++nonzero_count; }
    }
    msg_size += nonzero_count * sizeof(Vector);
  }

  void *msg_buf = CkAllocBuffer(msg,msg_size);
  char *msg_cur = (char *)msg_buf;

  CmiMemcpy((void*)msg_cur,(void*)(&(msg->node)),sizeof(msg->node));
  msg_cur += sizeof(msg->node);
  CmiMemcpy((void*)msg_cur,(void*)(&(msg->patch)),sizeof(msg->patch));
  msg_cur += sizeof(msg->patch);
  for ( j = 0; j < Results::maxNumForces; ++j ) {
    int array_size = msg->forceList[j].size();
    *(int *) msg_cur = array_size;
    msg_cur += sizeof(int);
    char *nonzero = msg_cur;
    msg_cur += array_size * sizeof(char);
    msg_cur = (char *)ALIGN_8 (msg_cur);
    Vector *farr = (Vector *)msg_cur;
    Force* f = msg->forceList[j].begin();

    for ( int i = 0; i < array_size; ++i ) {
      if ( f[i].x != 0. || f[i].y != 0. || f[i].z != 0. ) {
        nonzero[i] = 1;
	farr->x = f[i].x;
	farr->y = f[i].y;
	farr->z = f[i].z;
	farr ++;
      } else {
        nonzero[i] = 0;
      }
    }
    msg_cur = (char *) farr;	  
  }

  delete msg;
  return msg_buf;
}

ProxyResultMsg* ProxyResultMsg::unpack(void *ptr) {
  void *vmsg = CkAllocBuffer(ptr,sizeof(ProxyResultMsg));
  ProxyResultMsg *msg = new (vmsg) ProxyResultMsg;
  char *msg_cur = (char*)ptr;

  CmiMemcpy((void*)(&(msg->node)),(void*)msg_cur,sizeof(msg->node));
  msg_cur += sizeof(msg->node);
  CmiMemcpy((void*)(&(msg->patch)),(void*)msg_cur,sizeof(msg->patch));
  msg_cur += sizeof(msg->patch);
  int j;
  for ( j = 0; j < Results::maxNumForces; ++j ) {
    int array_size = *(int *) msg_cur;
    msg_cur += sizeof(array_size);
    msg->forceList[j].resize(array_size);
    char *nonzero = msg_cur;
    msg_cur += array_size * sizeof(char);    
    msg_cur = (char *)ALIGN_8 (msg_cur);
    Vector* farr = (Vector *) msg_cur;
    Force* f = msg->forceList[j].begin();
    for ( int i = 0; i < array_size; ++i ) {
      if ( nonzero[i] ) {
	f[i].x = farr->x;
	f[i].y = farr->y;
	f[i].z = farr->z;
	farr++;
      } else {
        f[i].x = 0.;  f[i].y = 0.;  f[i].z = 0.;
      }
    }    
    msg_cur = (char *) farr;
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
    msg_size = ALIGN_8 (msg_size);

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
  CmiMemcpy((void*)msg_cur,(void*)(&nodeSize), sizeof(nodeSize));
  msg_cur += sizeof(nodeSize);
  for (int i=0; i<nodeSize; i++) {
    CmiMemcpy((void*)msg_cur,(void*)(&msg->nodes[i]), sizeof(NodeID));
    msg_cur += sizeof(NodeID);
  }
  CmiMemcpy((void*)msg_cur,(void*)(&(msg->patch)),sizeof(msg->patch));
  msg_cur += sizeof(msg->patch);
  for ( j = 0; j < Results::maxNumForces; ++j ) {
    int array_size = msg->forceList[j].size();
    CmiMemcpy((void*)msg_cur,(void*)(&array_size),sizeof(array_size));
    msg_cur += sizeof(array_size);
    char *nonzero = msg_cur;
    msg_cur += array_size * sizeof(char);
    msg_cur = (char *)ALIGN_8 (msg_cur);
    Vector *farr = (Vector *) msg_cur; 
    Force* f = msg->forceList[j].begin();

    for ( int i = 0; i < array_size; ++i ) {
      if ( f[i].x != 0. || f[i].y != 0. || f[i].z != 0. ) {
        nonzero[i] = 1;
	farr->x  =  f[i].x;
        farr->y  =  f[i].y;
        farr->z  =  f[i].z;

        farr ++;
      } else {
        nonzero[i] = 0;
      }
    }
    msg_cur = (char *) farr;
  }

  delete msg;
  return msg_buf;
}

ProxyCombinedResultMsg* ProxyCombinedResultMsg::unpack(void *ptr) {
  void *vmsg = CkAllocBuffer(ptr,sizeof(ProxyCombinedResultMsg));
  ProxyCombinedResultMsg *msg = new (vmsg) ProxyCombinedResultMsg;
  char *msg_cur = (char*)ptr;

  int nodeSize;
  CmiMemcpy((void*)(&nodeSize),(void*)msg_cur,sizeof(nodeSize));
  msg_cur += sizeof(nodeSize);
  for (int i=0; i<nodeSize; i++) {
    msg->nodes.add(*(int *)msg_cur);
    msg_cur += sizeof(NodeID);
  }
  CmiMemcpy((void*)(&(msg->patch)),(void*)msg_cur,sizeof(msg->patch));
  msg_cur += sizeof(msg->patch);
  int j;
  for ( j = 0; j < Results::maxNumForces; ++j ) {
    int array_size;
    CmiMemcpy((void*)(&array_size),(void*)msg_cur,sizeof(array_size));
    msg_cur += sizeof(array_size);
    msg->forceList[j].resize(array_size);
    char *nonzero = msg_cur;
    msg_cur += array_size * sizeof(char);
    msg_cur = (char *)ALIGN_8 (msg_cur);
    Vector* farr = (Vector *) msg_cur;
    Force* f = msg->forceList[j].begin();

    for ( int i = 0; i < array_size; ++i ) {
      if ( nonzero[i] ) {
	f[i].x = farr->x;
	f[i].y = farr->y;
	f[i].z = farr->z;
	farr++;
      } else {
        f[i].x = 0.;  f[i].y = 0.;  f[i].z = 0.;
      }
    }
    msg_cur = (char *) farr;
  }

  CkFreeMsg(ptr);
  return msg;
}

// class static
int ProxyMgr::nodecount = 0;

ProxyMgr::ProxyMgr() { 
  if (CpvAccess(ProxyMgr_instance)) {
    NAMD_bug("Tried to create ProxyMgr twice.");
  }
  CkpvAccess(ProxyMgr_instance) = this;
}

ProxyMgr::~ProxyMgr() { 
  removeProxies();
  CpvAccess(ProxyMgr_instance) = NULL;
}


void ProxyMgr::setSendSpanning() {
  proxySendSpanning = 1;
}

int ProxyMgr::getSendSpanning() {
  return proxySendSpanning;
}

void ProxyMgr::setRecvSpanning() {
  proxySendSpanning = 1;
}

int ProxyMgr::getRecvSpanning() {
  return proxySendSpanning;
}

ProxyTree &ProxyMgr::getPtree() {
  return ptree;
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
      //fprintf(stderr, "Proxy Deleted Patch %d Proc %d", pi->patchID, CkMyPe());
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
  PatchIDList basepids;
  patchMap->basePatchIDList(myNode,basepids);
  for ( i = 0; i < basepids.size(); ++i )
  {
    if ( patchMap->node(basepids[i]) != myNode ) {
	patchFlag[basepids[i]] = NeedProxy;
    }
    int numNeighbors = patchMap->upstreamNeighbors(basepids[i],neighbors);
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
    // iout << iINFO << "Removing unused proxy " << pid << " on " << iPE << ".\n" << endi;
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
#if CHARM_VERSION > 050402
  cp[node].recvRegisterProxy(msg);
#else
  cp.recvRegisterProxy(msg,node);
#endif
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
#if CHARM_VERSION > 050402
  cp[node].recvUnregisterProxy(msg);
#else
  cp.recvUnregisterProxy(msg,node);
#endif
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
  if (!CkMyPe()) iout << iINFO << "Building spanning tree ... send: " << proxySendSpanning << " recv: " << proxyRecvSpanning << "\n" << endi;
  PatchMap::Object()->homePatchIDList(pids);
  for (int i=0; i<pids.size(); i++) {
    HomePatch *home = PatchMap::Object()->homePatch(pids[i]);
    if (home == NULL) CkPrintf("ERROR: homepatch NULL\n");
    home->buildSpanningTree();
  }
}

void 
ProxyMgr::buildProxySpanningTree2()
{
  PatchIDList pids;
  PatchMap::Object()->homePatchIDList(pids);
  for (int i=0; i<pids.size(); i++) {
    HomePatch *home = PatchMap::Object()->homePatch(pids[i]);
    if (home == NULL) CkPrintf("ERROR: homepatch NULL\n");
    home->sendProxies();
  }
}

void 
ProxyMgr::sendProxies(int pid, int *list, int n)
{
  CProxy_ProxyMgr cp(CpvAccess(BOCclass_group).proxyMgr);
  cp[0].recvProxies(pid, list, n);
}

// only on PE 0
void 
ProxyMgr::recvProxies(int pid, int *list, int n)
{
  int nPatches = PatchMap::Object()->numPatches();
  if (ptree.proxylist == NULL)
    ptree.proxylist = new NodeIDList[nPatches];
  ptree.proxylist[pid].resize(n);
  for (int i=0; i<n; i++)
    ptree.proxylist[pid][i] = list[i];
  ptree.proxyMsgCount ++;
  if (ptree.proxyMsgCount == nPatches) {
    ptree.proxyMsgCount = 0;
    // building and sending of trees is done in two steps now
    // so that the building step can be shifted to the load balancer
    buildSpanningTree0();
    sendSpanningTrees();
  }
}

#define MAX_INTERNODE 1

//
// XXX static and global variables are unsafe for shared memory builds.
// The global and static vars should be eliminated.  
// Unfortunately, the routines that use these below are actually 
// in use in NAMD.
//
extern double *cpuloads;
static int *procidx = NULL;
static double averageLoad = 0.0;

static int compLoad(const void *a, const void *b)
{
  int i1 = *(int *)a;
  int i2 = *(int *)b;
  double d1 = cpuloads[i1];
  double d2 = cpuloads[i2];
  if (d1 < d2) 
    return 1;
  else if (d1 == d2) 
    return 0;
  else 
    return -1;
  // sort from high to low
}

static void processCpuLoad()
{
  int i;
  if (!procidx) {
    procidx = new int[CkNumPes()];
  }
  for (i=0; i<CkNumPes(); i++) procidx[i] = i;
  qsort(procidx, CkNumPes(), sizeof(int), compLoad);

  double averageLoad = 0.0;
  for (i=0; i<CkNumPes(); i++) averageLoad += cpuloads[i];
  averageLoad /= CkNumPes();
//  iout << "buildSpanningTree1: no intermediate node on " << procidx[0] << " " << procidx[1] << endi;

}

static int noInterNode(int p)
{
  int exclude = 0;
  if(CkNumPes()<1025)
    exclude = 5;
  else if(CkNumPes()<4097)
    exclude = 10;
  else if(CkNumPes()<8193)
    exclude = 40;
  else if(CkNumPes()<16385)
    exclude = 40;
  else
    exclude = 80;
  for (int i=0; i<exclude; i++) if (procidx[i] == p) return 1;
//  if (cpuloads[p] > averageLoad) return 1;
  return 0;
}

// only on PE 0
void 
ProxyMgr::buildSpanningTree0()
{
  int i;

  processCpuLoad();

  int *numPatchesOnNode = new int[CkNumPes()];
  int numNodesWithPatches = 0;
  for (i=0; i<CkNumPes(); i++) numPatchesOnNode[i] = 0;
  int numPatches = PatchMap::Object()->numPatches();
  for (i=0; i<numPatches; i++) {
    int node = PatchMap::Object()->node(i);
    numPatchesOnNode[node]++;
    if (numPatchesOnNode[node] == 1)
      numNodesWithPatches ++;
  }
  int patchNodesLast =
    ( numNodesWithPatches < ( 0.7 * CkNumPes() ) );
  int *ntrees = new int[CkNumPes()];
  for (i=0; i<CkNumPes(); i++) ntrees[i] = 0;
  if (ptree.trees == NULL) ptree.trees = new NodeIDList[numPatches];
  for (int pid=0; pid<numPatches; pid++) 
  {
    int numProxies = ptree.proxylist[pid].size();
    if (numProxies == 0) {
      CkPrintf ("This is sheer evil!\n\n");
      //ProxyMgr::Object()->sendSpanningTreeToHomePatch(pid, NULL, 0);
      return;
    }
    NodeIDList &tree = ptree.trees[pid];   // spanning tree
    NodeIDList oldtree = tree;
    tree.resize(numProxies+1);
    tree.setall(-1);
    tree[0] = PatchMap::Object()->node(pid);
    int s=1, e=numProxies;
    int nNonPatch = 0;
    int treesize = 1;
    int pp;

    // keep tree persistent for non-intermediate nodes
    for (pp=0; pp<numProxies; pp++) {
      int p = ptree.proxylist[pid][pp];
      int oldindex = oldtree.find(p);
      if (oldindex != -1 && oldindex <= numProxies) {
        int isIntermediate = (oldindex*proxySpanDim+1 <= numProxies);
        if (!isIntermediate) {
          tree[oldindex] = p;
        }
        else if (ntrees[p] < MAX_INTERNODE) {
          tree[oldindex] = p;
          ntrees[p] ++;
        }
      }
    }

    for (pp=0; pp<numProxies; pp++) {
      int p = ptree.proxylist[pid][pp];              // processor number
      if (tree.find(p) != -1) continue;        // already used
      treesize++;
      if (patchNodesLast && numPatchesOnNode[p] ) {
        while (tree[e] != -1) { e--; if (e==-1) e = numProxies; }
        tree[e] = p;
        int isIntermediate = (e*proxySpanDim+1 <= numProxies);
        if (isIntermediate) ntrees[p]++;
      }
      else {
        while (tree[s] != -1) { s++; if (s==numProxies+1) s = 1; }
        int isIntermediate = (s*proxySpanDim+1 <= numProxies);
        if (isIntermediate && (ntrees[p] >= MAX_INTERNODE || noInterNode(p))) {   // TOO MANY INTERMEDIATE TREES
        //if (isIntermediate && ntrees[p] >= MAX_INTERNODE)    // TOO MANY INTERMEDIATE TREES
          while (tree[e] != -1) { e--; if (e==-1) e = numProxies; }
          tree[e] = p;
          isIntermediate = (e*proxySpanDim+1 <= numProxies);
          if (isIntermediate) ntrees[p]++;
        }
        else {
          tree[s] = p;
          nNonPatch++;
          if (isIntermediate) ntrees[p]++;
        }
      }
    }
    // send homepatch's proxy tree
    if(ptree.sizes)
      ptree.sizes[pid] = treesize;
    //ProxyMgr::Object()->sendSpanningTreeToHomePatch(pid, &tree[0], treesize);
  }
  /*for (i=0; i<CkNumPes(); i++) {
    if (ntrees[i] > MAX_INTERNODE) iout << "Processor " << i << "has (guess) " << ntrees[i] << " intermediate nodes." << endi;
  }*/
  delete [] ntrees;
  delete [] numPatchesOnNode;
}

void ProxyMgr::sendSpanningTrees()
{
  int numPatches = PatchMap::Object()->numPatches();
  for (int pid=0; pid<numPatches; pid++) {
    int numProxies = ptree.proxylist[pid].size();
    if (numProxies == 0)
      ProxyMgr::Object()->sendSpanningTreeToHomePatch(pid, NULL, 0);
    else {
      /*if(ptree.sizes)
        ProxyMgr::Object()->sendSpanningTreeToHomePatch(pid, ptree.trees[pid].begin(), ptree.sizes[pid]);
      else*/
        ProxyMgr::Object()->sendSpanningTreeToHomePatch(pid, ptree.trees[pid].begin(), ptree.trees[pid].size());
    }
  }
}

void ProxyMgr::sendSpanningTreeToHomePatch(int pid, int *tree, int n)
{
  CProxy_ProxyMgr cp(thisgroup);
  cp[PatchMap::Object()->node(pid)].recvSpanningTreeOnHomePatch(pid, tree, n);
}

void ProxyMgr::recvSpanningTreeOnHomePatch(int pid, int *tree, int n)
{
  HomePatch *p = PatchMap::Object()->homePatch(pid);
  p->recvSpanningTree(tree, n);
}

void 
ProxyMgr::sendSpanningTree(ProxySpanningTreeMsg *msg) {
  CProxy_ProxyMgr cp(CpvAccess(BOCclass_group).proxyMgr);
#if CHARM_VERSION > 050402
  cp[msg->tree[0]].recvSpanningTree(msg);
#else
  cp.recvSpanningTree(msg, msg->tree[0]);
#endif
}

void 
ProxyMgr::recvSpanningTree(ProxySpanningTreeMsg *msg) {
  int size = msg->tree.size();
  int child[proxySpanDim];
  int nChild = 0;
  int i;
  ProxyPatch *proxy = (ProxyPatch *) PatchMap::Object()->patch(msg->patch);
  for (i=0; i<proxySpanDim; i++) {
    if (size > i+1) { child[i] = msg->tree[i+1]; nChild++; }
  }
  if (!PatchMap::Object()->homePatch(msg->patch)) {
    proxy->setSpanningTree(msg->node, child, nChild);
  }

  // build subtree and pass down
  if (nChild == 0) return;

  nodecount ++;
  //if (nodecount > MAX_INTERNODE) 
  //  iout << "Processor " << CkMyPe() << "has (actual) " << nodecount << " intermediate nodes." << endi;

//CkPrintf("[%d] %d:(%d) %d %d %d %d %d\n", CkMyPe(), msg->patch, size, msg->tree[0], msg->tree[1], msg->tree[2], msg->tree[3], msg->tree[4]);
  NodeIDList *tree = new NodeIDList[proxySpanDim];
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
    level *=proxySpanDim;
  }

  ProxyMgr *proxyMgr = ProxyMgr::Object();
  for (i=0; i<proxySpanDim; i++) {
    if (tree[i].size()) {
      ProxySpanningTreeMsg *cmsg = new ProxySpanningTreeMsg;
      cmsg->patch = msg->patch;
      cmsg->node = CkMyPe();
      cmsg->tree = tree[i];
      proxyMgr->sendSpanningTree(cmsg);
    }
  }

  delete [] tree;
  delete msg;
}

void
ProxyMgr::sendResults(ProxyResultMsg *msg) {
  CProxy_ProxyMgr cp(CpvAccess(BOCclass_group).proxyMgr);
  NodeID node = PatchMap::Object()->node(msg->patch);
#if CHARM_VERSION > 050402
  cp[node].recvResults(msg);
#else
  cp.recvResults(msg, node);
#endif
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
    int destPe = patch->getSpanningTreeParent();
    if(destPe != CkMyPe()) {
#if CHARM_VERSION > 050402
      cp[destPe].recvImmediateResults(cMsg);
#else
      cp.recvImmediateResults(cMsg, destPe);
#endif
    }
    else 
      cp[destPe].recvResults(cMsg);
  }
}

void
ProxyMgr::recvResults(ProxyCombinedResultMsg *msg) {
  HomePatch *home = PatchMap::Object()->homePatch(msg->patch);
  if (home) {
    //printf("Home got a message\n");
    home->receiveResults(msg); // delete done in HomePatch::receiveResults()
  }
  else {
    NAMD_bug("ProxyMgr should receive result message on home processor");
  }
}

void ProxyMgr::recvImmediateResults(ProxyCombinedResultMsg *msg) {
  HomePatch *home = PatchMap::Object()->homePatch(msg->patch);
  if (home) {
    CProxy_ProxyMgr cp(CpvAccess(BOCclass_group).proxyMgr);
#if CHARM_VERSION > 050402
    cp[CkMyPe()].recvResults(msg);
#else
    cp.recvResults(msg, CkMyPe());
#endif
  }
  else {
    ProxyPatch *patch = (ProxyPatch *)PatchMap::Object()->patch(msg->patch);
    ProxyCombinedResultMsg *cMsg = patch->depositCombinedResultMsg(msg);
    if (cMsg) {
      CProxy_ProxyMgr cp(CpvAccess(BOCclass_group).proxyMgr);
#if CHARM_VERSION > 050402
      cp[patch->getSpanningTreeParent()].recvImmediateResults(cMsg);
#else
      cp.recvImmediateResults(cMsg, patch->getSpanningTreeParent());
#endif
    }
  }
}

void
ProxyMgr::sendProxyData(ProxyDataMsg *msg, int pcnt, int *pids) {
  CProxy_ProxyMgr cp(CpvAccess(BOCclass_group).proxyMgr);
  cp.recvImmediateProxyData(msg,pcnt,pids);
}

void 
ProxyMgr::recvProxyData(ProxyDataMsg *msg) {
  ProxyPatch *proxy = (ProxyPatch *) PatchMap::Object()->patch(msg->patch);
  proxy->receiveData(msg); // deleted in ProxyPatch::receiveAtoms()
}

void
ProxyMgr::recvImmediateProxyData(ProxyDataMsg *msg) {
  ProxyPatch *proxy = (ProxyPatch *) PatchMap::Object()->patch(msg->patch);
  if (proxySendSpanning == 1) {
    // copy the message and send to spanning children
    int pids[proxySpanDim];
    int npid = proxy->getSpanningTreeChild(pids);
    if (npid) {
      ProxyDataMsg *newmsg = new(PRIORITY_SIZE) ProxyDataMsg;
      CkSetQueueing(newmsg, CK_QUEUEING_IFIFO);
      *((int*) CkPriorityPtr(newmsg)) = *((int*) CkPriorityPtr(msg));
      newmsg->patch = msg->patch;
      newmsg->flags = msg->flags;
      newmsg->positionList = msg->positionList;
      newmsg->avgPositionList = msg->avgPositionList;
      //ProxyMgr::Object()->sendProxyData(newmsg,npid,pids);

      //At the second level of the tree immediate messages are not needed
      CProxy_ProxyMgr cp(CpvAccess(BOCclass_group).proxyMgr);
      cp.recvProxyData(newmsg,npid,pids);
    }
  }
  /* send to self via EP method to preserve priority */
  CProxy_ProxyMgr cp(CpvAccess(BOCclass_group).proxyMgr);
  cp[CkMyPe()].recvProxyData(msg);
}

void
ProxyMgr::sendProxyAll(ProxyAllMsg *msg, int pcnt, int *pids) {
  CProxy_ProxyMgr cp(CpvAccess(BOCclass_group).proxyMgr);
  cp.recvImmediateProxyAll(msg,pcnt,pids);
}

void 
ProxyMgr::recvProxyAll(ProxyAllMsg *msg) {
  ProxyPatch *proxy = (ProxyPatch *) PatchMap::Object()->patch(msg->patch);
  proxy->receiveAll(msg); // deleted in ProxyPatch::receiveAtoms()
}

void
ProxyMgr::recvImmediateProxyAll(ProxyAllMsg *msg) {
  ProxyPatch *proxy = (ProxyPatch *) PatchMap::Object()->patch(msg->patch);
  if (proxySendSpanning == 1) {
    // copy the message and send to spanning children
    int pids[proxySpanDim];
    int npid = proxy->getSpanningTreeChild(pids);
    if (npid) {
      ProxyAllMsg *newmsg = new(PRIORITY_SIZE) ProxyAllMsg;
      CkSetQueueing(newmsg, CK_QUEUEING_IFIFO);
      *((int*) CkPriorityPtr(newmsg)) = *((int*) CkPriorityPtr(msg));
      newmsg->patch = msg->patch;
      newmsg->flags = msg->flags;
      newmsg->positionList = msg->positionList;
      newmsg->avgPositionList = msg->avgPositionList;
#ifdef MEM_OPT_VERSION
      newmsg->extInfoList = msg->extInfoList;
#endif
      ProxyMgr::Object()->sendProxyAll(newmsg,npid,pids);
    }
  }
  /* send to self via EP method to preserve priority */
  CProxy_ProxyMgr cp(CpvAccess(BOCclass_group).proxyMgr);
  cp[CkMyPe()].recvProxyAll(msg);
}

#include "ProxyMgr.def.h"

