/***************************************************************************/
/*                                                                         */
/*              (C) Copyright 1995 The Board of Trustees of the            */
/*                          University of Illinois                         */
/*                           All Rights Reserved                           */
/*									   */
/***************************************************************************/
 
/***************************************************************************
 * DESCRIPTION:
 *	ComputeSphericalBC provides spherical boundary conditions.  See the
 * Programmer's Guide for details on the potentials provided.
 *
 ***************************************************************************/

#ifndef COMPUTESPHERICALBC_H
#define COMPUTESPHERICALBC_H

#include "ComputePatch.h"
#include "ReductionMgr.h"

class ComputeSphericalBC : public ComputePatch
{
private:
	BigReal r1;			//  Radius of first sphere
	BigReal r1_2;			//  Radius of first sphere squared
	BigReal k1;			//  First force constant
	BigReal r2;			//  Radius of second sphere (-1 if inactive)
	BigReal r2_2;			//  Raidus of second sphere squared
	BigReal k2;			//  Second force constant
	int exp1;			//  Exponent for first boundary condition
	int exp2;			//  Exponent for second boundary condition
	Bool twoForces;			//  Are there two potentials or just one
	BigReal energy;			//  Energy computed for the current timestep
	Vector center;			//  Center of spheres

public:
	ComputeSphericalBC(ComputeID c, PatchID pid); 	//  Constructor
	virtual ~ComputeSphericalBC();			//  Destructor

	virtual void doForce(Position* p, Results* r, AtomProperties* a);

	SubmitReduction *reduction;

};

#endif







