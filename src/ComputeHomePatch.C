/**
***  Copyright (c) 1995, 1996, 1997, 1998, 1999, 2000 by
***  The Board of Trustees of the University of Illinois.
***  All rights reserved.
**/

/*
   Compute object which deals with a single patch.
*/

#include "charm++.h"
#include "WorkDistrib.decl.h"
#include "Node.h"
#include "ComputeHomePatch.h"
#include "PatchMap.inl"
#include "HomePatch.h"

#define MIN_DEBUG_LEVEL 4
//#define DEBUGM
#include "Debug.h"

ComputeHomePatch::ComputeHomePatch(ComputeID c, PatchID p) : Compute(c) {
    setNumPatches(1);
    patchID = p;
    patch = NULL;
    homePatch = NULL;
    positionBox = NULL;
    forceBox = NULL;
    atomBox = NULL;
}

ComputeHomePatch::~ComputeHomePatch() {
  DebugM(4, "~ComputeHomePatch("<<cid<<") numAtoms("<<patchID<<") = " 
    << numAtoms << "\n");
    if (positionBox != NULL) {
      PatchMap::Object()->patch(patchID)->unregisterPositionPickup(cid,
	 &positionBox);
    }
    if (forceBox != NULL) {
      PatchMap::Object()->patch(patchID)->unregisterForceDeposit(cid,
		&forceBox);
    }
    if (atomBox != NULL) {
      PatchMap::Object()->patch(patchID)->unregisterAtomPickup(cid,
		&atomBox);
    }
}

void ComputeHomePatch::initialize() {
    // How can we tell if BoxOwner has packed up and left?  Need a mechanism
    // to handle this or do we assume the Boxes have been dumped?

	if (positionBox == NULL) { // We have yet to get boxes
	    if (!(patch = PatchMap::Object()->patch(patchID))) {
	      NAMD_bug("ComputeHomePatch used with unknown patch.");
	    }
            if (!(homePatch = PatchMap::Object()->homePatch(patchID))) {
	      NAMD_bug("ComputeHomePatch used with proxy.");
	    }
	    DebugM(3, "initialize(" << cid <<")  patchid = "<<patch->getPatchID()<<"\n");
	    positionBox = patch->registerPositionPickup(cid);
	    forceBox = patch->registerForceDeposit(cid);
	    atomBox = patch->registerAtomPickup(cid);
	}
	numAtoms = patch->getNumAtoms();

  DebugM(3, "initialize("<<cid<<") numAtoms("<<patchID<<") = " 
    << numAtoms  << " patchAddr=" << patch << "\n");
    Compute::initialize();

    int myNode = CkMyPe();
    if ( PatchMap::Object()->node(patchID) != myNode )
    {
      basePriority = patchID % 64;
    }
    else
    {
      basePriority = 2 * 64 + (patchID % 64);
    }
}

void ComputeHomePatch::atomUpdate() {
    // How can we tell if BoxOwner has packed up and left?  Need a mechanism
    // to handle this or do we assume the Boxes have been dumped?
    numAtoms = patch->getNumAtoms();
}

void ComputeHomePatch::doWork() {
  Position* p;
  Results* r;
  AtomProperties* a;
  Transform* t = homePatch->getTransformList().begin();
  int numData;

  DebugM(3,patchID << ": doWork() called.\n");

  // Open up positionBox, forceBox, and atomBox
  p = positionBox->open(&numData);
  if (numData != numAtoms) {
    NAMD_bug("doWork has opened a position box with wrong # atoms.");
  }
  r = forceBox->open();
  a = atomBox->open();

  // Pass pointers to doForce
  doForce(p,r,a,t);

  // Close up boxes
  positionBox->close(&p);
  forceBox->close(&r);
  atomBox->close(&a);

  DebugM(2,patchID << ": doWork() completed.\n");
}

int ComputeHomePatch::sequence(void)
{
  return patch->flags.step;
}

