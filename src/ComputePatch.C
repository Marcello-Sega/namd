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

#include "WorkDistrib.top.h"
#include "Node.h"
#include "ComputePatch.h"

#define MIN_DEBUG_LEVEL 4
#define DEBUGM
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

void ComputePatch::mapReady() {
    // How can we tell if BoxOwner has packed up and left?  Need a mechanism
    // to handle this or do we assume the Boxes have been dumped?

	if (positionBox == NULL) { // We have yet to get boxes
	    patch = PatchMap::Object()->patch(patchID);
	    positionBox = patch->registerPositionPickup(cid);
	    forceBox = patch->registerForceDeposit(cid);
	    atomBox = patch->registerAtomPickup(cid);
	}
	numAtoms = patch->getNumAtoms();

    Compute::mapReady();
}

void ComputePatch::doForce(Position* p,
                               Force* f,
                               AtomProperties* a)
{
    CPrintf("ComputePatchPair::doForce() - Dummy eval was sent\n");
    CPrintf(" %d patch 1 atoms\n", numAtoms );
}

void ComputePatch::doWork() {
  Position* p;
  Force* f;
  AtomProperties* a;
  int i;

  DebugM(3,patchID << ": doWork() called.\n");

  // Open up positionBox, forceBox, and atomBox
      p = positionBox->open();
      f = forceBox->open();
      a = atomBox->open();

  // Pass pointers to doForce
  doForce(p,f,a);

  // Close up boxes
      DebugM(1,patchID << ": closing positionBox.\n");
      positionBox->close(&p);
      DebugM(1,patchID << ": closing forceBox.\n");
      forceBox->close(&f);
      DebugM(1,patchID << ": closing atomBox.\n");
      atomBox->close(&a);

  DebugM(2,patchID << ": doWork() completed.\n");

}


/***************************************************************************
 * RCS INFORMATION:
 *
 *	$RCSfile: ComputePatch.C,v $
 *	$Author: jim $	$Locker:  $		$State: Exp $
 *	$Revision: 1.6 $	$Date: 1996/12/01 02:39:58 $
 *
 ***************************************************************************
 * REVISION HISTORY:
 *
 * $Log: ComputePatch.C,v $
 * Revision 1.6  1996/12/01 02:39:58  jim
 * added debugging
 *
 * Revision 1.5  1996/10/31 22:05:55  jim
 * first incarnation as ComputePatch
 *
 * Revision 1.4  1996/10/30 01:16:32  jim
 * added AtomProperties structure in Patch plus boxes, passing, etc.
 *
 * Revision 1.3  1996/10/30 00:16:16  jim
 * Removed PositionArray usage.
 *
 * Revision 1.2  1996/10/29 23:53:58  jim
 * cleaned up, now only compile blocks are PatchMap, Patch, Compute.
 *
 * Revision 1.1  1996/10/29 22:43:35  ari
 * Initial revision
 *
 * Revision 1.2  1996/10/22 19:12:16  ari
 * *** empty log message ***
 *
 * Revision 1.1  1996/10/16 08:22:39  ari
 * Initial revision
 *
 * Revision 1.1  1996/08/19 22:07:49  ari
 * Initial revision
 *
 * Revision 1.4  1996/07/16 01:54:12  ari
 * *** empty log message ***
 *
 * Revision 1.3  96/07/16  01:10:26  01:10:26  ari (Aritomo Shinozaki)
 * Fixed comments, added methods
 * 
 * Revision 1.2  1996/06/25 21:10:48  gursoy
 * *** empty log message ***
 *
 * Revision 1.1  1996/06/24 14:12:26  gursoy
 * Initial revision
 *
 ***************************************************************************/

