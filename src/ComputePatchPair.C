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

#include "WorkDistrib.top.h"
#include "Node.h"
#include "ComputePatchPair.h"
#include "Patch.h"

//#define DEBUGM
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
	    patch[i] = PatchMap::Object()->patch(patchID[i]);
	    positionBox[i] = patch[i]->registerPositionPickup(cid,trans[i]);
	    forceBox[i] = patch[i]->registerForceDeposit(cid);
	    atomBox[i] = patch[i]->registerAtomPickup(cid);
	}
	numAtoms[i] = patch[i]->getNumAtoms();
    }

    Compute::initialize();
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
    CPrintf("ComputePatchPair::doForce() - Dummy eval was sent\n");
    CPrintf(" %d patch 1 atoms and %d patch 2 atoms\n", numAtoms[0], numAtoms[1] );
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

  // Open up positionBox, forceBox, and atomBox
  for (i=0; i<2; i++) {
      p[i] = positionBox[i]->open();
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


/***************************************************************************
 * RCS INFORMATION:
 *
 *	$RCSfile: ComputePatchPair.C,v $
 *	$Author: jim $	$Locker:  $		$State: Exp $
 *	$Revision: 1.1007 $	$Date: 1997/03/18 18:08:56 $
 *
 ***************************************************************************
 * REVISION HISTORY:
 *
 * $Log: ComputePatchPair.C,v $
 * Revision 1.1007  1997/03/18 18:08:56  jim
 * Revamped collection system to ensure ordering and eliminate
 * unnecessary collections.  Also reduced make dependencies.
 *
 * Revision 1.1006  1997/03/13 06:37:09  jim
 * Multiple time-stepping implemented, still needs proper splitting functions.
 *
 * Revision 1.1005  1997/03/12 22:06:38  jim
 * First step towards multiple force returns and multiple time stepping.
 *
 * Revision 1.1004  1997/03/06 22:06:00  ari
 * Removed Compute.ci
 * Comments added - more code cleaning
 *
 * Revision 1.1003  1997/02/13 23:17:17  ari
 * Fixed a final bug in AtomMigration - numatoms in ComputePatchPair.C not
 * set correctly in atomUpdate()
 *
 * Revision 1.1002  1997/02/11 18:51:45  ari
 * Modified with #ifdef DPMTA to safely eliminate DPMTA codes
 * fixed non-buffering of migration msgs
 * Migration works on multiple processors
 *
 * Revision 1.1001  1997/02/07 07:50:23  jim
 * Removed debug messages.
 *
 * Revision 1.1000  1997/02/06 15:58:17  ari
 * Resetting CVS to merge branches back into the main trunk.
 * We will stick to main trunk development as suggested by CVS manual.
 * We will set up tags to track fixed points of development/release
 * as suggested by CVS manual - all praise the CVS manual.
 *
 * Revision 1.779  1997/02/06 15:53:09  ari
 * Updating Revision Line, getting rid of branches
 *
 * Revision 1.778.2.2  1997/02/06 02:35:21  jim
 * Implemented periodic boundary conditions - may not work with
 * atom migration yet, but doesn't seem to alter calculation,
 * appears to work correctly when turned on.
 * NamdState chdir's to same directory as config file in argument.
 *
 * Revision 1.778.2.1  1997/02/05 22:18:11  ari
 * Added migration code - Currently the framework is
 * there with compiling code.  This version does
 * crash shortly after migration is complete.
 * Migration appears to complete, but Patches do
 * not appear to be left in a correct state.
 *
 * Revision 1.778  1997/01/28 00:30:29  ari
 * internal release uplevel to 1.778
 *
 * Revision 1.777  1997/01/17 19:36:03  ari
 * Internal CVS leveling release.  Start development code work
 * at 1.777.1.1.
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

