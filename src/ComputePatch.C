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
#include "Patch.h"

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

void ComputePatch::initialize() {
    // How can we tell if BoxOwner has packed up and left?  Need a mechanism
    // to handle this or do we assume the Boxes have been dumped?

	if (positionBox == NULL) { // We have yet to get boxes
	    patch = PatchMap::Object()->patch(patchID);
	    positionBox = patch->registerPositionPickup(cid);
	    forceBox = patch->registerForceDeposit(cid);
	    atomBox = patch->registerAtomPickup(cid);
	}
	numAtoms = patch->getNumAtoms();

    Compute::initialize();
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

  DebugM(3,patchID << ": doWork() called.\n");

  // Open up positionBox, forceBox, and atomBox
  p = positionBox->open();
  r = forceBox->open();
  a = atomBox->open();

  // Pass pointers to doForce
  doForce(p,r,a);

  // Close up boxes
  positionBox->close(&p);
  forceBox->close(&r);
  atomBox->close(&a);

  DebugM(2,patchID << ": doWork() completed.\n");
}


/***************************************************************************
 * RCS INFORMATION:
 *
 *	$RCSfile: ComputePatch.C,v $
 *	$Author: jim $	$Locker:  $		$State: Exp $
 *	$Revision: 1.1004 $	$Date: 1997/03/18 18:08:54 $
 *
 ***************************************************************************
 * REVISION HISTORY:
 *
 * $Log: ComputePatch.C,v $
 * Revision 1.1004  1997/03/18 18:08:54  jim
 * Revamped collection system to ensure ordering and eliminate
 * unnecessary collections.  Also reduced make dependencies.
 *
 * Revision 1.1003  1997/03/13 06:37:06  jim
 * Multiple time-stepping implemented, still needs proper splitting functions.
 *
 * Revision 1.1002  1997/03/12 22:06:37  jim
 * First step towards multiple force returns and multiple time stepping.
 *
 * Revision 1.1001  1997/02/11 18:51:44  ari
 * Modified with #ifdef DPMTA to safely eliminate DPMTA codes
 * fixed non-buffering of migration msgs
 * Migration works on multiple processors
 *
 * Revision 1.1000  1997/02/06 15:58:15  ari
 * Resetting CVS to merge branches back into the main trunk.
 * We will stick to main trunk development as suggested by CVS manual.
 * We will set up tags to track fixed points of development/release
 * as suggested by CVS manual - all praise the CVS manual.
 *
 * Revision 1.779  1997/02/06 15:53:07  ari
 * Updating Revision Line, getting rid of branches
 *
 * Revision 1.778.2.1  1997/02/05 22:18:10  ari
 * Added migration code - Currently the framework is
 * there with compiling code.  This version does
 * crash shortly after migration is complete.
 * Migration appears to complete, but Patches do
 * not appear to be left in a correct state.
 *
 * Revision 1.778  1997/01/28 00:30:27  ari
 * internal release uplevel to 1.778
 *
 * Revision 1.777  1997/01/17 19:36:01  ari
 * Internal CVS leveling release.  Start development code work
 * at 1.777.1.1.
 *
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

