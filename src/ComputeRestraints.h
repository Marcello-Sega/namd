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
	//****** BEGIN selective restraints (X,Y,Z) changes 
	Bool consSelectOn;      // Selection of Cartesian components active?
	Bool consSelectX, consSelectY,
	     consSelectZ;       // which components are active?
	//****** END selective restraints (X,Y,Z) changes 
	//****** BEGIN moving constraints changes 
	Bool consMoveOn;        //  Are the moving constraints on?
        Vector moveVel;         // velocity of the constraint movement (A/timestep).
	//****** END moving constraints changes 
	//****** BEGIN rotating constraints changes 
	// rotating constraints. 
	// Ref. pos. of all the atoms that are constrained will rotate
	Bool consRotOn;         // Are the rotating constraints on?
	Vector rotAxis;         // Axis of rotation
        Vector rotPivot;        // Pivot point of rotation
        BigReal rotVel;         // Rotation velocity (deg/timestep);
	//****** END rotating constraints changes 

public:
	ComputeRestraints(ComputeID c, PatchID pid); 	//  Constructor
	virtual ~ComputeRestraints();			//  Destructor

	virtual void doForce(Position* p, Results* r, AtomProperties* a);

	SubmitReduction *reduction;

};

#endif

/***************************************************************************
 * RCS INFORMATION:
 *
 *	$RCSfile: ComputeRestraints.h,v $
 *	$Author: jim $	$Locker:  $		$State: Exp $
 *	$Revision: 1.5 $	$Date: 1999/06/17 15:46:13 $
 *
 ***************************************************************************
 * REVISION HISTORY:
 *
 * $Log: ComputeRestraints.h,v $
 * Revision 1.5  1999/06/17 15:46:13  jim
 * Completely rewrote reduction system to eliminate need for sequence numbers.
 *
 * Revision 1.4  1999/01/08 23:24:34  ferenc
 * added selective position restraints for specific Cartesian components
 *
 * Revision 1.3  1998/10/01 00:31:31  sergei
 * added rotating restraints feature;
 * changed the moving restraints from only moving one atom to moving all
 * atoms that are restrained. One-atom pulling is available in SMD feature.
 *
 * Revision 1.2  1997/08/18 20:16:06  sergei
 * added moving restraint capability with input from config file
 * one atom only
 *
 * Revision 1.1  1997/04/22 04:26:01  jim
 * Added atomic restraints (harmonic constraints) via ComputeRestraints class.
 *
 *
 ***************************************************************************/






