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

#include "main.h"
#include "ckdefs.h"
#include "chare.h"
#include "c++interface.h"

#include "WorkDistrib.top.h"

#include "NamdTypes.h"
#include "Templates/Box.h"
#include "Templates/OwnerBox.h"

#include "Node.h"
#include "Compute.h"

ComputePatchPair::ComputePatchPair(ComputeID c, PatchID p[2]) : Compute(c) {
  setNumPatches(2);
  for (int i=0; i<2; i++) {
      patchID[i] = p[i];
      patch[i] = NULL;
      positionBox[i] = NULL;
      forceBox[i] = NULL;
  }
}

ComputePatchPair::~ComputePatchPair() {
  for (int i=0; i<2; i++) {
    if (positionBox[i] != NULL) {
      PatchMap::global->patch(pid[i])->unregisterPositionPickup(cid,
	 &positionBox[i]);
    }
    if (positionBox[i] != NULL) {
      PatchMap::global->patch(pid[i])->unregisterForceDeposit(cid,
		&forceBox[i]);
    }
  }

}

void ComputePatchPair::mapReady() {
    // How can we tell if BoxOwner has packed up and left?  Need a mechanism
    // to handle this or do we assume the Boxes have been dumped?

    for (int i=0; i<2; i++) {
	if (positionBox[i] == NULL) { // We have yet to get boxes
	    patch[i] = PatchMap::global->patch(pid[i]);
	    positionBox[i] = patch[i]->registerPositionPickup(cid);
	    forceBox[i] = patch[i]->registerForceDeposit(cid);
	}
	numAtoms[i] = patch[i]->getNumAtoms();
    }
}

void ComputePatchPair::depositAllForces() {
  for (i=0; i<2; i++) {
      positionBox[i]->close();
      forceBox[i]->close();
  }
}

void ComputePatchPair::doForce(Position p[2][], Force f[2][]) {
    CPrintf("ComputePatchPair::doForce() - Dummy eval was sent\n");
    CPrintf(" %d patch 1 atoms and %d patch 2 atoms\n", numAtoms[0], numAtoms[1] );
}

void ComputePatchPair::doWork() {
  Position *p[2];
  Force *f[2];

  // Open up positionBox and forceBox
  for (int i=0; i<2; i++) {
      p[i] = positionBox[i]->open();
      f[i] = forceBox[i]->open();
  }

  // Pass pointers to doForce
  doForce( p, f, numAtoms );

  // Close up boxes
  // Open up positionBox and forceBox
  for (int i=0; i<2; i++) {
      p[i] = positionBox[i]->close();
      f[i] = forceBox[i]->close();
  }
}


/***************************************************************************
 * RCS INFORMATION:
 *
 *	$RCSfile: ComputePatch.C,v $
 *	$Author: ari $	$Locker:  $		$State: Exp $
 *	$Revision: 1.1 $	$Date: 1996/10/29 22:43:35 $
 *
 ***************************************************************************
 * REVISION HISTORY:
 *
 * $Log: ComputePatch.C,v $
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

