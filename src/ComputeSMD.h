//-*-c++-*-
/***************************************************************************/
/*                                                                         */
/*              (C) Copyright 1995 The Board of Trustees of the            */
/*                          University of Illinois                         */
/*                           All Rights Reserved                           */
/*									   */
/***************************************************************************/
 
/***************************************************************************
 * DESCRIPTION:
 *	ComputeSMD provides moving atomic restraints (harmonic constraints).  
 * See the Programmer's Guide for details on the potentials provided.
 *
 ***************************************************************************/

#ifndef COMPUTESMD_H
#define COMPUTESMD_H

#include "ComputePatch.h"
#include "ReductionMgr.h"

class ComputeSMD : public ComputePatch
{
private:
	int consExp;		//  Exponent for energy function from SimParameters
        BigReal k;              //  Restraint force constant
	int moveAtom;           //  Index of the atom to move
        BigReal moveVel;         // velocity of the restraint movement 
                                // (A/timestep).
        int outputFreq;         // output frequency
        Bool chDirOn;           // is changing direction on?
        Bool chForceOn;         // is changing force on?

public:
	ComputeSMD(ComputeID c, PatchID pid); 	//  Constructor
	virtual ~ComputeSMD();			//  Destructor

	virtual void doForce(Position* p, Results* r, AtomProperties* a);

	ReductionMgr *reduction;

};

#endif







