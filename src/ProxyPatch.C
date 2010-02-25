/**
***  Copyright (c) 1995, 1996, 1997, 1998, 1999, 2000 by
***  The Board of Trustees of the University of Illinois.
***  All rights reserved.
**/

#include "InfoStream.h"
#include "Lattice.h"
#include "main.decl.h"
#include "main.h"
#include "ProxyPatch.h"
#include "ProxyMgr.decl.h"
#include "ProxyMgr.h"
#include "AtomMap.h"
#include "PatchMap.h"
#include "Priorities.h"

#define MIN_DEBUG_LEVEL 2
//#define  DEBUGM
#include "Debug.h"


ProxyPatch::ProxyPatch(PatchID pd) : 
  Patch(pd), proxyMsgBufferStatus(PROXYMSGNOTBUFFERED), 
  curProxyMsg(NULL), prevProxyMsg(NULL)
{
  DebugM(4, "ProxyPatch(" << pd << ") at " << this << "\n");
  ProxyMgr::Object()->registerProxy(patchID);
  numAtoms = -1;
  parent = -1;

#ifdef NODEAWARE_PROXY_SPANNINGTREE
  /*numChild = 0;
  children = NULL;*/
#else
  nChild = 0;
  child = new int[proxySpanDim];
#endif

#if CMK_PERSISTENT_COMM
  localphs = 0;
  localphs = CmiCreatePersistent(PatchMap::Object()->node(patchID), 300000);
#endif

  // DMK - Atom Separation (water vs. non-water)
  #if NAMD_SeparateWaters != 0
    numWaterAtoms = -1;
  #endif
  
  #if defined(NODEAWARE_PROXY_SPANNINGTREE) && defined(USE_NODEPATCHMGR)
    depositLock = CmiCreateLock();
  #endif
}

ProxyPatch::~ProxyPatch()
{
  DebugM(4, "ProxyPatch(" << patchID << ") deleted at " << this << "\n");
  ProxyMgr::Object()->unregisterProxy(patchID);

  // ProxyPatch may be freed because of load balancing if the compute object
  // it corresponds to no longer exist on this specific processor.
  CmiAssert(prevProxyMsg!=NULL);
  if(prevProxyMsg!=NULL) {
// #ifdef REMOVE_PROXYDATAMSG_EXTRACOPY
//       AtomMap::Object()->unregisterIDs(patchID,positionPtrBegin, positionPtrEnd);
// #else
      AtomMap::Object()->unregisterIDs(patchID,pExt.begin(),pExt.end());
// #endif      
      delete prevProxyMsg;
      prevProxyMsg = NULL;
  }


#ifdef NODEAWARE_PROXY_SPANNINGTREE
  delete [] children;
  #ifdef USE_NODEPATCHMGR
  delete [] nodeChildren;  
  #endif
#else
  delete [] child;
#endif

  p.resize(0);
  pExt.resize(0);

#if CMK_PERSISTENT_COMM
  CmiDestoryPersistent(localphs);
  localphs = 0;
#endif
}

void ProxyPatch::boxClosed(int box)
{
  if ( box == 1 ) {	
    // Note: delay the deletion of proxyDataMsg (of the 
    // current step) until the next step. This is done 
    // for the sake of atom migration (ProxyDataMsg) 
    // as the ProxyPatch has to  unregister the atoms 
    // of the previous step in the AtomMap data structure. 
    sendResults();
  }
  if ( ! --boxesOpen ) {
    DebugM(2,patchID << ": " << "Checking message buffer.\n");    
    
    if(proxyMsgBufferStatus == PROXYALLMSGBUFFERED) {
          CmiAssert(curProxyMsg != NULL);
          DebugM(3,"Patch " << patchID << " processing buffered proxy ALL data.\n");
          receiveAll(curProxyMsg);          
    }else if(proxyMsgBufferStatus == PROXYDATAMSGBUFFERED) {
          CmiAssert(curProxyMsg != NULL);
          DebugM(3,"Patch " << patchID << " processing buffered proxy data.\n");
          receiveData(curProxyMsg);
    }
  } else {
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

  //delete the ProxyDataMsg of the previous step
  delete prevProxyMsg;
  prevProxyMsg = NULL;

  if ( boxesOpen )
  {
      proxyMsgBufferStatus = PROXYDATAMSGBUFFERED;
    // store message in queue (only need one element, though)
    curProxyMsg = msg;
    return;
  }

  //Reuse position arrays inside proxyDataMsg --Chao Mei
  curProxyMsg = msg;
  prevProxyMsg = curProxyMsg;
  flags = msg->flags;

#ifdef REMOVE_PROXYDATAMSG_EXTRACOPY
  //We could set them to 0 for the sake of easy debugging
  //if there are something wrong in the "reuse position arrays" code
  //--Chao Mei
  //p.resize(0);
  //p_avg.resize(0);  
  positionPtrBegin = msg->positionList;
  positionPtrEnd = msg->positionList + msg->plLen;
#else
  p.resize(msg->plLen);
  memcpy(p.begin(), msg->positionList, sizeof(CompAtom)*(msg->plLen));
#endif
  
  avgPositionPtrBegin = msg->avgPositionList;
  avgPositionPtrEnd = msg->avgPositionList + msg->avgPlLen;
  
  // BEGIN LA
  velocityPtrBegin = msg->velocityList;
  velocityPtrEnd = msg->velocityList + msg->vlLen;
  // END LA

  if ( numAtoms == -1 ) { // for new proxies since receiveAtoms is not called
      //numAtoms = p.size();
      numAtoms = msg->plLen;

      //Retrieve the CompAtomExt list
      CmiAssert(msg->plExtLen!=0);
      pExt.resize(msg->plExtLen);
      memcpy(pExt.begin(), msg->positionExtList, sizeof(CompAtomExt)*(msg->plExtLen));


    // DMK - Atom Separation (water vs. non-water)
    #if NAMD_SeparateWaters != 0
      numWaterAtoms = msg->numWaterAtoms;
    #endif

    positionsReady(1);
  } else {
    positionsReady(0);
  }
}

void ProxyPatch::receiveAll(ProxyDataMsg *msg)
{
  DebugM(3, "receiveAll(" << patchID << ")\n");

  if ( boxesOpen )
  {
    proxyMsgBufferStatus = PROXYALLMSGBUFFERED;    
    curProxyMsg = msg;
    return;
  }  

  //The prevProxyMsg has to be deleted after this if-statement because
  // positionPtrBegin points to the space inside the prevProxyMsg
  if(prevProxyMsg!=NULL) {
// #ifdef REMOVE_PROXYDATAMSG_EXTRACOPY
//       AtomMap::Object()->unregisterIDs(patchID,positionPtrBegin,positionPtrEnd);
// #else
      AtomMap::Object()->unregisterIDs(patchID, pExt.begin(), pExt.end());
// #endif
  }
  //Now delete the ProxyDataMsg of the previous step
  delete prevProxyMsg;
  curProxyMsg = msg;
  prevProxyMsg = curProxyMsg;

  flags = msg->flags;

#ifdef REMOVE_PROXYDATAMSG_EXTRACOPY
  //We could set them to 0 for the sake of easy debugging
  //if there are something wrong in the "reuse position arrays" code
  //--Chao Mei
  //p.resize(0);
  //p_avg.resize(0);  
  positionPtrBegin = msg->positionList;
  positionPtrEnd = msg->positionList + msg->plLen;
#else
  p.resize(msg->plLen);
  memcpy(p.begin(), msg->positionList, sizeof(CompAtom)*(msg->plLen));
#endif

  numAtoms = msg->plLen;
  //numAtoms = p.size();
  
  avgPositionPtrBegin = msg->avgPositionList;
  avgPositionPtrEnd = msg->avgPositionList + msg->avgPlLen;
  
  // BEGIN LA
  velocityPtrBegin = msg->velocityList;
  velocityPtrEnd = msg->velocityList + msg->vlLen;
  // END LA

  //We cannot reuse the CompAtomExt list inside the msg because
  //the information is needed at every step. In the current implementation
  //scheme, the ProxyDataMsg msg will be deleted for every step.
  //In order to keep this information, we have to do the extra copy. But
  //this overhead is amortized among the steps that atoms don't migrate
  // --Chao Mei
  pExt.resize(msg->plExtLen);
  memcpy(pExt.begin(), msg->positionExtList, sizeof(CompAtomExt)*(msg->plExtLen));

  // DMK - Atom Separation (water vs. non-water)
  #if NAMD_SeparateWaters != 0
    numWaterAtoms = msg->numWaterAtoms;
  #endif

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
#ifdef REMOVE_PROXYRESULTMSG_EXTRACOPY
    ProxyResultVarsizeMsg *msg = ProxyResultVarsizeMsg::getANewMsg(CkMyPe(), patchID, PRIORITY_SIZE, f); 
#else
    ProxyResultMsg *msg = new (PRIORITY_SIZE) ProxyResultMsg;    
    msg->node = CkMyPe();
    msg->patch = patchID;
    for ( i = 0; i < Results::maxNumForces; ++i ) 
      msg->forceList[i] = f[i];
#endif
    SET_PRIORITY(msg,flags.sequence,PROXY_RESULTS_PRIORITY + PATCH_PRIORITY(patchID));
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

#ifdef NODEAWARE_PROXY_SPANNINGTREE
void ProxyPatch::setSpanningTree(int p, int *c, int n) { 
  parent=p; numChild = n; nWait = 0;
  delete [] children;
  if(n==0) {
      children = NULL;
      return;
  }
  children = new int[n];
  for (int i=0; i<n; i++) children[i] = c[i];

  #if defined(PROCTRACE_DEBUG) && defined(NAST_DEBUG)
    DebugFileTrace *dft = DebugFileTrace::Object();
    dft->openTrace();
    dft->writeTrace("ProxyPatch[%d] has %d children: ", patchID, numChild);
    for(int i=0; i<numChild; i++)
        dft->writeTrace("%d ", children[i]);
    dft->writeTrace("\n");
    dft->closeTrace();
  #endif
//CkPrintf("setSpanningTree: [%d:%d] %d %d:%d %d\n", CkMyPe(), patchID, parent, nChild, child[0], child[1]);
}

int ProxyPatch::getSpanningTreeChild(int *c) { 
  for (int i=0; i<numChild; i++) c[i] = children[i];
  return numChild;
}

#ifdef USE_NODEPATCHMGR
void ProxyPatch::setSTNodeChildren(int numNids, int *nids){
    numNodeChild = numNids;
    delete [] nodeChildren;
    if(numNids==0) {
        nodeChildren = NULL;
        return;
    }
    nodeChildren = new int[numNids];
    for(int i=0; i<numNids; i++) nodeChildren[i] = nids[i]; 
}
#endif

#else //branch for not defined NODEAWARE_PROXY_SPANNINGTREE
void ProxyPatch::setSpanningTree(int p, int *c, int n) { 
  parent=p; nChild = n; nWait = 0;
  for (int i=0; i<n; i++) child[i] = c[i];
//CkPrintf("setSpanningTree: [%d:%d] %d %d:%d %d\n", CkMyPe(), patchID, parent, nChild, child[0], child[1]);
}

int ProxyPatch::getSpanningTreeChild(int *c) { 
  for (int i=0; i<nChild; i++) c[i] = child[i];
  return nChild;
}
#endif

ProxyCombinedResultMsg *ProxyPatch::depositCombinedResultMsg(ProxyCombinedResultMsg *msg) {
#if defined(NODEAWARE_PROXY_SPANNINGTREE) && defined(USE_NODEPATCHMGR)
  CmiLock(depositLock);
#endif
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
#ifdef NODEAWARE_PROXY_SPANNINGTREE
  if(nWait == numChild+1) {
#else
  if (nWait == nChild + 1) {
#endif
    nWait = 0;
#if defined(NODEAWARE_PROXY_SPANNINGTREE) && defined(USE_NODEPATCHMGR)
    CmiUnlock(depositLock);
#endif
    
    return msgCBuffer;
  }

#if defined(NODEAWARE_PROXY_SPANNINGTREE) && defined(USE_NODEPATCHMGR)
  CmiUnlock(depositLock);
#endif

  return NULL;
}

