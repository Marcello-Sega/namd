/***************************************************************************/
/*                                                                         */
/*              (C) Copyright 1996 The Board of Trustees of the            */
/*                          University of Illinois                         */
/*                           All Rights Reserved                           */
/*									   */
/***************************************************************************/

/***************************************************************************
 * DESCRIPTION: ComputeNonbondedPair.C
 *
 ***************************************************************************/

#include "ComputeNonbondedPair.h"

void ComputeNonbondedPair::doForce(Position* p[2],
                               Force* f[2],
                               AtomProperties* a[2])
{
    CPrintf("ComputeNonbondedPair::doForce() - Dummy eval was sent\n");
    CPrintf(" %d patch 1 atoms and %d patch 2 atoms\n", numAtoms[0], numAtoms[1] );
}

/***************************************************************************
 * RCS INFORMATION:
 *
 *	$RCSfile: ComputeNonbondedPair.C,v $
 *	$Author: jim $	$Locker:  $		$State: Exp $
 *	$Revision: 1.4 $	$Date: 1996/10/31 21:43:29 $
 *
 ***************************************************************************
 * REVISION HISTORY:
 *
 * $Log: ComputeNonbondedPair.C,v $
 * Revision 1.4  1996/10/31 21:43:29  jim
 * First incarnation as ...Pair
 *
 * Revision 1.3  1996/10/30 01:16:32  jim
 * added AtomProperties structure in Patch plus boxes, passing, etc.
 *
 * Revision 1.2  1996/10/30 00:16:16  jim
 * Removed PositionArray usage.
 *
 * Revision 1.1  1996/10/29 23:55:54  jim
 * Initial revision
 *
 *
 ***************************************************************************/

