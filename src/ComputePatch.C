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
      PatchMap::Object()->patch(patchID[i])->unregisterPositionPickup(cid,
	 &positionBox[i]);
    }
    if (positionBox[i] != NULL) {
      PatchMap::Object()->patch(patchID[i])->unregisterForceDeposit(cid,
		&forceBox[i]);
    }
  }

}

void ComputePatchPair::mapReady() {
    // How can we tell if BoxOwner has packed up and left?  Need a mechanism
    // to handle this or do we assume the Boxes have been dumped?

    for (int i=0; i<2; i++) {
	if (positionBox[i] == NULL) { // We have yet to get boxes
	    patch[i] = PatchMap::Object()->patch(patchID[i]);
	    positionBox[i] = patch[i]->registerPositionPickup(cid);
	    forceBox[i] = patch[i]->registerForceDeposit(cid);
	}
	numAtoms[i] = patch[i]->getNumAtoms();
    }

    Compute::mapReady();
}

/*
void ComputePatchPair::depositAllForces() {
  for (int i=0; i<2; i++) {
      positionBox[i]->close();
      forceBox[i]->close();
  }
}
*/

void ComputePatchPair::doForce(PositionArray p[2], ForceArray f[2]) {
    CPrintf("ComputePatchPair::doForce() - Dummy eval was sent\n");
    CPrintf(" %d patch 1 atoms and %d patch 2 atoms\n", numAtoms[0], numAtoms[1] );
}

void ComputePatchPair::doWork() {
  Position *p[2];
  Force *f[2];
  int i;

  // Open up positionBox and forceBox
  for (i=0; i<2; i++) {
      p[i] = positionBox[i]->open();
      f[i] = forceBox[i]->open();
  }

  // Pass pointers to doForce
  doForce(p,f);

  // Close up boxes
  // Open up positionBox and forceBox
  for (i=0; i<2; i++) {
      positionBox[i]->close(&p[i]);
      forceBox[i]->close(&f[i]);
  }
}


/***************************************************************************
 * RCS INFORMATION:
 *
 *	$RCSfile: ComputePatch.C,v $
 *	$Author: jim $	$Locker:  $		$State: Exp $
 *	$Revision: 1.2 $	$Date: 1996/10/29 23:53:58 $
 *
 ***************************************************************************
 * REVISION HISTORY:
 *
 * $Log: ComputePatch.C,v $
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

