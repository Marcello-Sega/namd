/***************************************************************************/
/*                                                                         */
/*              (C) Copyright 1996 The Board of Trustees of the            */
/*                          University of Illinois                         */
/*                           All Rights Reserved                           */
/*									   */
/***************************************************************************/

/***************************************************************************
 * DESCRIPTION: ComputeNonbonded.C
 *
 ***************************************************************************/

#include "ComputeNonbonded.h"

void ComputeNonbonded::doForce(PositionArray p[2], ForceArray f[2]) {
    CPrintf("ComputeNonbonded::doForce() - Dummy eval was sent\n");
    CPrintf(" %d patch 1 atoms and %d patch 2 atoms\n", numAtoms[0], numAtoms[1] );
}

/***************************************************************************
 * RCS INFORMATION:
 *
 *	$RCSfile: ComputeNonbondedPair.C,v $
 *	$Author: jim $	$Locker:  $		$State: Exp $
 *	$Revision: 1.1 $	$Date: 1996/10/29 23:55:54 $
 *
 ***************************************************************************
 * REVISION HISTORY:
 *
 * $Log: ComputeNonbondedPair.C,v $
 * Revision 1.1  1996/10/29 23:55:54  jim
 * Initial revision
 *
 *
 ***************************************************************************/

