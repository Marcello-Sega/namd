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

