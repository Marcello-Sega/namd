/***************************************************************************/
/*                                                                         */
/*              (C) Copyright 1995 The Board of Trustees of the            */
/*                          University of Illinois                         */
/*                           All Rights Reserved                           */
/*									   */
/***************************************************************************/
 
/***************************************************************************
 * DESCRIPTION:
 *	ComputeRestraints provides atomic restraints (harmonic constraints).  See the
 * Programmer's Guide for details on the potentials provided.
 *
 ***************************************************************************/

#ifndef COMPUTERESTRAINTS_H
#define COMPUTERESTRAINTS_H

#include "ComputePatch.h"
#include "ReductionMgr.h"

class ComputeRestraints : public ComputePatch
{
private:
	int consExp;		//  Exponent for energy function from SimParameters
	Bool consMoveOn;        //  Are the cmoving constraints on?
	int moveAtom;           //  Index of the atom to move
        Vector moveVel;         // velocity of the constraint movement (A/timestep).

public:
	ComputeRestraints(ComputeID c, PatchID pid); 	//  Constructor
	virtual ~ComputeRestraints();			//  Destructor

	virtual void doForce(Position* p, Results* r, AtomProperties* a);

	ReductionMgr *reduction;

};

#endif

/***************************************************************************
 * RCS INFORMATION:
 *
 *	$RCSfile: ComputeRestraints.h,v $
 *	$Author: sergei $	$Locker:  $		$State: Exp $
 *	$Revision: 1.2 $	$Date: 1997/08/18 20:16:06 $
 *
 ***************************************************************************
 * REVISION HISTORY:
 *
 * $Log: ComputeRestraints.h,v $
 * Revision 1.2  1997/08/18 20:16:06  sergei
 * added moving restraint capability with input from config file
 * one atom only
 *
 * Revision 1.1  1997/04/22 04:26:01  jim
 * Added atomic restraints (harmonic constraints) via ComputeRestraints class.
 *
 *
 ***************************************************************************/






