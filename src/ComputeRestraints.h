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
 *	$Author: jim $	$Locker:  $		$State: Exp $
 *	$Revision: 1.1 $	$Date: 1997/04/22 04:26:01 $
 *
 ***************************************************************************
 * REVISION HISTORY:
 *
 * $Log: ComputeRestraints.h,v $
 * Revision 1.1  1997/04/22 04:26:01  jim
 * Added atomic restraints (harmonic constraints) via ComputeRestraints class.
 *
 *
 ***************************************************************************/






