/***************************************************************************/
/*                                                                         */
/*              (C) Copyright 1996 The Board of Trustees of the            */
/*                          University of Illinois                         */
/*                           All Rights Reserved                           */
/*									   */
/***************************************************************************/

/***************************************************************************
 * DESCRIPTION: Compute object which deals with a single patch.
 *
 ***************************************************************************/

#include "charm++.h"
#include "WorkDistrib.decl.h"
#include "Node.h"
#include "ComputePatch.h"
#include "PatchMap.inl"
#include "Patch.h"

#define MIN_DEBUG_LEVEL 4
//#define DEBUGM
#include "Debug.h"

ComputePatch::ComputePatch(ComputeID c, PatchID p) : Compute(c) {
    setNumPatches(1);
    patchID = p;
    patch = NULL;
    positionBox = NULL;
    forceBox = NULL;
    atomBox = NULL;
}

ComputePatch::~ComputePatch() {
  DebugM(4, "~ComputePatch("<<cid<<") numAtoms("<<patchID<<") = " 
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

void ComputePatch::initialize() {
    // How can we tell if BoxOwner has packed up and left?  Need a mechanism
    // to handle this or do we assume the Boxes have been dumped?

	if (positionBox == NULL) { // We have yet to get boxes
	    if (!(patch = PatchMap::Object()->patch(patchID))) {
	      iout << iERRORF << iPE << "invalid patch! during initialize\n" 
		   << endi;
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

void ComputePatch::atomUpdate() {
    // How can we tell if BoxOwner has packed up and left?  Need a mechanism
    // to handle this or do we assume the Boxes have been dumped?
    numAtoms = patch->getNumAtoms();
}

void ComputePatch::doForce(Position* p,
                               Results* r,
                               AtomProperties* a)
{
    DebugM(1, "ComputePatchPair::doForce() - Dummy eval was sent\n");
    if (p && r && a) {
      p[0] = Position(0.0,0.0,0.0);
    }
}

void ComputePatch::doWork() {
  Position* p;
  Results* r;
  AtomProperties* a;
  int numData;

  DebugM(3,patchID << ": doWork() called.\n");

  // Open up positionBox, forceBox, and atomBox
  p = positionBox->open(&numData);
  if (numData != numAtoms) {
    iout << iPE << iERRORF 
      << "Interesting, doWork has opened a position box with wrong # atoms ("
      <<numData<<" vs " << numAtoms << "\n" 
      << endi;
  }
  r = forceBox->open();
  a = atomBox->open();

  /*
  if (!p || !r || !a) {
    iout << iPE << iERRORF
     << "Data Pointer is NULL! on open on patchID(" << patchID
     << ")\n" << endi;
  }
  */

  // Pass pointers to doForce
  doForce(p,r,a);

  // Close up boxes
  positionBox->close(&p);
  forceBox->close(&r);
  atomBox->close(&a);

  DebugM(2,patchID << ": doWork() completed.\n");
}

int ComputePatch::sequence(void)
{
  return patch->flags.step;
}

