
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
#include "Sync.h"

#include "PatchMap.h"
#include "Compute.h"
#include "ComputeMap.h"

#define MIN_DEBUG_LEVEL 3
//#define DEBUGM
#include "Debug.h"

#include <stdio.h>
#include <converse.h>

#include "InfoStream.h"

int useSync = 1;

Sync::Sync()
{
    if (CpvAccess(Sync_instance) == NULL) {
        CpvAccess(Sync_instance) = this;
    } else {
	iout << iFILE << iERROR << iPE
	  << "Sync instanced twice on same processor!" << endi;
	CkExit();
    }
    counter = 0;
}

void Sync::PatchReady(void)
{
   counter ++;
//   CkPrintf("SYNC[%d]: PATCHREADY:%d %d \n", CkMyPe(), counter, PatchMap::Object()->numHomePatches());
   if (counter == PatchMap::Object()->numHomePatches()) 
   {
       counter = 0;
       // inform all computes to ready
       ComputeMap *computeMap = ComputeMap::Object();
       for (int i=0; i<computeMap->numComputes(); i++) {
	   if (computeMap->node(i) == CkMyPe()) {
	       computeMap->compute(i)->patchReady(0,0);
	   }
       }
   }
}


#include "Sync.def.h"
