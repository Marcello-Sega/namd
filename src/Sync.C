
/**
***  Copyright (c) 1995, 1996, 1997, 1998, 1999, 2000 by
***  The Board of Trustees of the University of Illinois.
***  All rights reserved.
**/

/*
   Toplevel routines for initializing a Node for a simulation
   one Node per Pe (processor element).
*/

#ifndef WIN32
#include <unistd.h>
#endif
#include <charm++.h>
#include "Sync.decl.h"

#include "Patch.h"
#include "PatchMap.h"
#include "ProxyMgr.h"
#include "Compute.h"
#include "ComputeMap.h"

#include "Sync.h"

#define MIN_DEBUG_LEVEL 3
//#define DEBUGM
#include "Debug.h"

#include <stdio.h>
#include <converse.h>

#include "InfoStream.h"

int useSync = 0;

Sync::Sync()
{
    if (CpvAccess(Sync_instance) == NULL) {
        CpvAccess(Sync_instance) = this;
    } else {
	iout << iFILE << iERROR << iPE
	  << "Sync instanced twice on same processor!" << endi;
	CkExit();
    }
    capacity = 3000;
    clist = new _clist[capacity];
    step = 0;
    counter = 0;
    cnum = 0;
    nPatcheReady = 0;
    numPatches = -1;
}

void Sync::openSync(void)
{
   if (useSync && PatchMap::Object()->numHomePatches() == 0) {
       CkPrintf("********* Local Sync is removed on node %d ******** \n", CkMyPe());
       useSync = 0;
   }
}    

void Sync::registerComp(PatchID pid, ComputeIDListIter cid, int doneMigration)
{
  int i;
  int slot = 0;
  for (slot = 0; slot < cnum; slot++)
     if (clist[slot].pid == -1) break;
  if (slot == cnum) {
    cnum++;
    // expand the list
    if (cnum == capacity) {
      capacity += 1000;
      struct _clist *tmp = new _clist[capacity];
      for (i=0; i<cnum; i++) tmp[i] = clist[i];
      delete [] clist;
      clist = tmp;
      CmiPrintf("Sync buffer overflow and expanded!\n");
    }
  }

  PatchMap *patchMap =  PatchMap::Object();
  clist[slot].cid = cid;
  clist[slot].pid = pid;
  clist[slot].doneMigration  = doneMigration;
  clist[slot].step = patchMap->patch(pid)->flags.step;

//  CkPrintf("REG[%d]: patch:%d step:%d-%d slot:%d\n", CkMyPe(), pid, patchMap->patch(pid)->flags.step, step, slot);

  triggerCompute();
}

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

  nPatcheReady= 0 ;
  for (int i= 0; i<cnum; i++) {
      PatchID pid = clist[i].pid;
      if (pid == -1) continue;
      if (clist[i].step == step) nPatcheReady++;
  }

  CkPrintf("SYNC[%d]: PATCHREADY:%d %d patches:%d %d\n", CkMyPe(), counter, PatchMap::Object()->numHomePatches(), nPatcheReady, numPatches);
  if (counter == patchMap->numHomePatches() && nPatcheReady == numPatches)
  {
       CkPrintf("TRIGGERED[%d]\n", CkMyPe());
       ComputeMap *computeMap = ComputeMap::Object();
       for (int i= 0; i<cnum; i++) {
         int &pid = clist[i].pid;
	 if (pid == -1) continue;
	 if (clist[i].step != step) continue;
         CkPrintf(" %d-%d-%d ", clist[i].pid, clist[i].step, patchMap->patch(pid)->flags.step);
         ComputeIDListIter cid = clist[i].cid;
         for(cid = cid.begin(); cid != cid.end(); cid++)
         {
	    computeMap->compute(*cid)->patchReady(pid,clist[i].doneMigration);
         }
	 pid = -1;
       }
       CkPrintf("\n");

       // reset counter
       counter = 0;
       numPatches = -1;
       step++;
  }
}


#include "Sync.def.h"
