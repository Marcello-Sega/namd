/***************************************************************************/
/*    (C) Copyright 1995,1996,1997 The Board of Trustees of the            */
/*                          University of Illinois                         */
/*                           All Rights Reserved                           */
/***************************************************************************/
/***************************************************************************
 * DESCRIPTION:
 *	ComputeCylindricalBC provides cylindrical boundary conditions.  See the
 * Programmer's Guide for details on the potentials provided.
 *
 ***************************************************************************/

#ifndef COMPUTECYLINDRICALBC_H
#define COMPUTECYLINDRICALBC_H

#include "ComputePatch.h"
#include "ReductionMgr.h"

class ComputeCylindricalBC : public ComputePatch
{
private:
	char axis;			//  'x', 'y', or 'z'
	BigReal r1;			//  Radius of first cylinder
	BigReal r1_2;			//  Radius of first cylinder squared
	BigReal l1;			//  Length of First cylinder
	BigReal l1_2;			//  Length of first cylinder, squared
	BigReal k1;			//  First force constant
	BigReal r2;			//  Radius of second cylinder (-1 if inactive)
	BigReal r2_2;			//  Raidus of second cylinder squared
	BigReal k2;			//  Second force constant
	BigReal l2;			//  Length of second cylinder
	BigReal l2_2;			//  Length of second cylinder, squared
	int exp1;			//  Exponent for first boundary condition
	int exp2;			//  Exponent for second boundary condition
	Bool twoForces;			//  Are there two potentials or just one
	Vector center;			//  Center of cylinder

public:
	ComputeCylindricalBC(ComputeID c, PatchID pid); 	//  Constructor
	virtual ~ComputeCylindricalBC();			//  Destructor

	virtual void doForce(Position* p, Results* r, AtomProperties* a);

	SubmitReduction *reduction;

};

#endif


/***************************************************************************
 * RCS INFORMATION:
 *
 *	$RCSfile $
 *	$Author $	$Locker:  $		$State: Exp $
 *	$Revision: 1.6 $	$Date: 1999/06/17 15:46:02 $
 *
 ***************************************************************************
 * REVISION HISTORY:
 *
 * $Log: ComputeCylindricalBC.h,v $
 * Revision 1.6  1999/06/17 15:46:02  jim
 * Completely rewrote reduction system to eliminate need for sequence numbers.
 *
 * Revision 1.5  1997/03/20 23:53:32  ari
 * Some changes for comments. Copyright date additions.
 * Hooks for base level update of Compute objects from ComputeMap
 * by ComputeMgr.  Useful for new compute migration functionality.
 *
 * Revision 1.4  1997/03/19 11:54:05  ari
 * Add Broadcast mechanism.
 * Fixed RCS Log entries on files that did not have Log entries.
 * Added some register variables to Molecule and ComputeNonbondedExcl.C
 *
 ***************************************************************************/
