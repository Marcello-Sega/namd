
/**
***  Copyright (c) 1995, 1996, 1997, 1998, 1999, 2000 by
***  The Board of Trustees of the University of Illinois.
***  All rights reserved.
**/

/*
    Sync will ensure that all homepatches finished updating before Computes starts and all proxies finished updating themselves.
*/

#if !defined(WIN32) || defined(__CYGWIN__)
#include <unistd.h>
#endif
#include <stdio.h>
#include <charm++.h>

#include "Patch.h"
#include "PatchMap.h"
#include "ProxyMgr.h"
#include "Compute.h"
#include "ComputeMap.h"

#include "Sync.h"

#define MIN_DEBUG_LEVEL 3
//#define DEBUGM
#include "Debug.h"

#include "InfoStream.h"

extern int proxySendSpanning, proxyRecvSpanning;

// make sure all HomePatches get their positions data and sendProxyData to 
// their proxies before computes get positionsReady.
int useSync = 1;

// useProxySync will make sure all proxies get updated before computes' 
// positionsReady triggered and real computation begins.
// when these two combined, it will make sure that homepatch get all its force 
// and positions data, and proxies receive its updated data before all 
// computes start.
int useProxySync = 0;

Sync::Sync()
{
    if (CpvAccess(Sync_instance) == NULL) {
        CpvAccess(Sync_instance) = this;
    } else {
	iout << iFILE << iERROR << iPE
	  << "Sync instanced twice on same processor!" << endi;
	CkExit();
    }
    capacity = 1000;
    clist = new _clist[capacity];
    step = 0;
    counter = 0;
    cnum = 0;
    nPatcheReady = 0;
    numPatches = -1;
}

Sync::~Sync()
{
  delete [] clist;
}

void Sync::openSync(void)
{
  if (useSync) {
    // if use proxy spanning tree, proxy sync is forced
    if (!useProxySync && (proxySendSpanning || proxyRecvSpanning)) useProxySync = 1;
    // no proxies on this node, no need to use proxy sync.
    if (useProxySync && ProxyMgr::Object()->numProxies() == 0) useProxySync = 0;
    // if no proxy sync and no home patch, then disable home patch sync as well
    if (!useProxySync && PatchMap::Object()->numHomePatches() == 0) useSync = 0;
  }
#if 0
// was
   if (useSync && PatchMap::Object()->numHomePatches() == 0) {
       // CkPrintf("********* Local Sync is removed on node %d ******** \n", CkMyPe());
     useSync = 0;
   }
#endif
}    

// called from Patch::positionsReady()
void Sync::registerComp(PatchID pid, ComputeIDListIter cid, int doneMigration)
{
  int slot = 0;
  for (; slot < cnum; slot++)
     if (clist[slot].pid == -1) break;
  if (slot == cnum) {
    cnum++;
    // table is full, expand the list
    if (cnum == capacity) {
      capacity += 1000;
      struct _clist *tmp = new _clist[capacity];
      memcpy(tmp, clist, cnum*sizeof(_clist));
      delete [] clist;
      clist = tmp;
      CmiPrintf("Info:: Sync buffer overflow and expanded!\n");
    }
  }

  clist[slot].cid = cid;
  clist[slot].pid = pid;
  clist[slot].doneMigration  = doneMigration;
  clist[slot].step = PatchMap::Object()->patch(pid)->flags.sequence;

//  CkPrintf("REG[%d]: patch:%d step:%d-%d slot:%d\n", CkMyPe(), pid, patchMap->patch(pid)->flags.sequence, step, slot);

  if (clist[slot].step == step) {
      nPatcheReady++;
      triggerCompute();
  }
}

// called from HomePatch::positionsReady()
void Sync::PatchReady(void)
{
  counter ++;
  triggerCompute();
}

void Sync::triggerCompute()
{
  PatchMap *patchMap = PatchMap::Object();

  if (numPatches == -1) 
    numPatches = ProxyMgr::Object()->numProxies() + patchMap->numHomePatches();

//  CkPrintf("SYNC[%d]: PATCHREADY:%d %d patches:%d %d\n", CkMyPe(), counter, PatchMap::Object()->numHomePatches(), nPatcheReady, numPatches);
  if (counter == patchMap->numHomePatches() && nPatcheReady == numPatches)
  {
//       CkPrintf("TRIGGERED[%d]\n", CkMyPe());
       ComputeMap *computeMap = ComputeMap::Object();
       nPatcheReady = 0;
       for (int i= 0; i<cnum; i++) {
         int &pid = clist[i].pid;
	 if (pid == -1) continue;
	 if (clist[i].step != step) {
            // count for next step
            if (clist[i].step == step + 1) nPatcheReady++;
            continue;
         }
	 //         CkPrintf(" %d-%d-%d ",
	 //	 clist[i].pid, clist[i].step,
	 //      patchMap->patch(pid)->flags.sequence);

         ComputeIDListIter cid = clist[i].cid;

         int compute_count = 0;
         for(cid = cid.begin(); cid != cid.end(); cid++) {
	   compute_count++;
	   computeMap->compute(*cid)->patchReady(pid,clist[i].doneMigration,step);
	 }
	 if (compute_count == 0 && patchMap->node(pid) != CkMyPe()) {
	   iout << iINFO << "PATCH_COUNT-Sync step " << step
		<< "]: Patch " << pid << " on PE " 
		<< CkMyPe() <<" home patch " 
		<< patchMap->node(pid) << " does not have any computes\n" 
		<< endi;
	 }
	 pid = -1;
       }
//       CkPrintf("\n");

       // reset counter
       counter = 0;
       numPatches = -1;
       step++;
  }
}


#include "Sync.def.h"
