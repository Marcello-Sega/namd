/***************************************************************************/
/*                                                                         */
/*              (C) Copyright 1996 The Board of Trustees of the            */
/*                          University of Illinois                         */
/*                           All Rights Reserved                           */
/*									   */
/***************************************************************************/

/***************************************************************************
 * DESCRIPTION:
 *
 ***************************************************************************/

#include "charm++.h"
#include "WorkDistrib.decl.h"
#include "Node.h"
#include "ComputePatchPair.h"
#include "PatchMap.inl"
#include "Patch.h"

//#define DEBUGM
#define MIN_DEBUG_LEVEL 4
#include "Debug.h"

ComputePatchPair::ComputePatchPair(ComputeID c, PatchID p[], int t[]) 
    : Compute(c) {

  setNumPatches(2);

  for (int i=0; i<2; i++) {
      patchID[i] = p[i];
      trans[i] = t[i];
      patch[i] = NULL;
      positionBox[i] = NULL;
      forceBox[i] = NULL;
      atomBox[i] = NULL;
  }
}

ComputePatchPair::~ComputePatchPair() {
  DebugM(4, "~ComputePatchPair("<<cid<<") numAtoms("<<patchID[0]<<") = " 
    << numAtoms[0] 
    << " numAtoms("<<patchID[1]<<") = " << numAtoms[1] << "\n" );
  DebugM(4, "~ComputePatchPair("<<cid<<") addr("<<patchID[0]<<") = " 
    << PatchMap::Object()->patch(patchID[0]) << " addr("<<patchID[1]<<") = "
    << PatchMap::Object()->patch(patchID[1]) << "\n");
  for (int i=0; i<2; i++) {
    if (positionBox[i] != NULL) {
      PatchMap::Object()->patch(patchID[i])->unregisterPositionPickup(cid,
	 &positionBox[i]);
    }
    if (forceBox[i] != NULL) {
      PatchMap::Object()->patch(patchID[i])->unregisterForceDeposit(cid,
		&forceBox[i]);
    }
    if (atomBox[i] != NULL) {
      PatchMap::Object()->patch(patchID[i])->unregisterAtomPickup(cid,
		&atomBox[i]);
    }
  }

}

void ComputePatchPair::initialize() {
    // How can we tell if BoxOwner has packed up and left?  Need a mechanism
    // to handle this or do we assume the Boxes have been dumped?

    for (int i=0; i<2; i++) {
	if (positionBox[i] == NULL) { // We have yet to get boxes
	    if (!(patch[i] = PatchMap::Object()->patch(patchID[i]))) {
	      iout << iPE << iERRORF << "invalid patch(" << patchID[i] 
		   << ")  pointer!\n" << endi;
	    }
	    positionBox[i] = patch[i]->registerPositionPickup(cid,trans[i]);
	    forceBox[i] = patch[i]->registerForceDeposit(cid);
	    atomBox[i] = patch[i]->registerAtomPickup(cid);
	}
	numAtoms[i] = patch[i]->getNumAtoms();
    }

  DebugM(4, "initialize("<<cid<<") numAtoms("<<patchID[0]<<") = " 
    << numAtoms[0] 
    << " numAtoms(" <<patchID[1]<<") = " << numAtoms[1] << "\n" );

    Compute::initialize();

    int myNode = CkMyPe();
    int p0 = patchID[0] % 64;
    int p1 = patchID[1] % 64;
    int patchPrio = ((p0<p1)?p0:p1);
    if ( PatchMap::Object()->node(patchID[0]) != myNode )
    {
      basePriority = patchPrio;
    }
    else if ( PatchMap::Object()->node(patchID[1]) != myNode )
    {
      basePriority = patchPrio;
    }
    else
    {
      basePriority = 2 * 64 + patchPrio;
    }

}

void ComputePatchPair::atomUpdate() {
    // How can we tell if BoxOwner has packed up and left?  Need a mechanism
    // to handle this or do we assume the Boxes have been dumped?

    // DebugM(4,"atomUpdate() - positionBox[0] is " << positionBox[0] << "\n");
    for (int i=0; i<2; i++) {
	numAtoms[i] = patch[i]->getNumAtoms();
    }

    // Compute::atomUpdate();
}



void ComputePatchPair::doForce(Position* p[2],
                               Results* r[2],
                               AtomProperties* a[2])
{
    CkPrintf("ComputePatchPair::doForce() - Dummy eval was sent\n");
    CkPrintf(" %d patch 1 atoms and %d patch 2 atoms\n", numAtoms[0], numAtoms[1] );
}

//---------------------------------------------------------------------
// Where the actual computation is invoked.  doForce is 
// overloaded with specific calculation
//---------------------------------------------------------------------
void ComputePatchPair::doWork() {
  Position* p[2];
  Results* r[2];
  AtomProperties* a[2];
  int i;
  int numData;

  // Open up positionBox, forceBox, and atomBox
  for (i=0; i<2; i++) {
      p[i] = positionBox[i]->open(&numData);
      if (numData != numAtoms[i]) {
	iout << iPE << iERRORF 
	  << "Interesting, doWork has opened a position box with wrong # atoms ("
	  <<numData<<" vs " << numAtoms << "\n" 
	  << endi;
      }
      r[i] = forceBox[i]->open();
      a[i] = atomBox[i]->open();
  }

  // Pass pointers to doForce
  doForce(p,r,a);

  // Close up boxes
  for (i=0; i<2; i++) {
      positionBox[i]->close(&p[i]);
      forceBox[i]->close(&r[i]);
      atomBox[i]->close(&a[i]);
  }
}

int ComputePatchPair::sequence(void)
{
  return patch[0]->flags.step;
}

