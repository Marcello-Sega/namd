/**
***  Copyright (c) 1995, 1996, 1997, 1998, 1999, 2000 by
***  The Board of Trustees of the University of Illinois.
***  All rights reserved.
**/

/*
   Top of Compute hierarchy.  
   enqueueWork() - delivers Compute object itself to queue up for doWork()
   doWork() - called by work queue
*/

#include "main.h"
#include "charm++.h"

#include "WorkDistrib.decl.h"
#include "WorkDistrib.h"

#include "NamdTypes.h"
#include "Box.h"
#include "OwnerBox.h"

#include "Node.h"
#include "Compute.h"

#include "Sync.h"

#define MIN_DEBUG_LEVEL 4
// #define DEBUGM
#include "Debug.h"

Node *Compute::node=0;
PatchMap *Compute::patchMap=0;

int Compute::totalComputes = 0;

Compute::Compute(ComputeID c) : basePriority(DEFPRIO), cid(c),
	localWorkMsg(new (sizeof(int)*8) LocalWorkMsg) { 
  totalComputes++;
  doAtomUpdate = false;
  computeType = ComputeMap::Object()->type(c);
}

Compute::~Compute() {
  delete localWorkMsg;
  totalComputes--;
}

void Compute::enqueueWork() {
  if (!this) { DebugM(4,"This Compute is NULL!!!\n"); }
  if ( ! noWork() )
  {
    WorkDistrib::messageEnqueueWork(this);  // should be in ComputeMgr?
  }
}

//---------------------------------------------------------------------
// Signal from patch or proxy that data is ready.
// When all Patches and Proxies needed by this Compute object
// have checked-in, we are ready to enqueueWork()
//---------------------------------------------------------------------
void Compute::patchReady(PatchID patchID, int doneMigration) { 
  if (doneMigration) { // If any patch has done migration - we must remap
    doAtomUpdate = true; 
  }

  if (numPatches <= 0) {
      DebugM(5,"Compute::patchReady("<<patchID<<")-call not valid!\n");
  } else {
//  gzheng
//    if (! --patchReadyCounter) {
    if ( --patchReadyCounter == (useSync?-1:0) ) {
      patchReadyCounter = numPatches;
      if (doAtomUpdate) {
	atomUpdate();
	doAtomUpdate = false;
      }
      enqueueWork();
    }
  }
}


int Compute::noWork() {
  return 0;
}

void Compute::doWork() {
    DebugM(5,"Default Compute::doWork() called.\n");
}

int Compute::sequence(void)
{
  return -1;
}

