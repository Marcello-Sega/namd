/***************************************************************************/
/*                                                                         */
/*              (C) Copyright 1995 The Board of Trustees of the            */
/*                          University of Illinois                         */
/*                           All Rights Reserved                           */
/*                                                                         */
/***************************************************************************/
       
#include "ComputeRestraints.h"
#include "Node.h"
#include "Molecule.h"
#include "SimParameters.h"
#include "Patch.h"
#include "NamdOneTools.h"

/************************************************************************/
/*									*/
/*			FUNCTION ComputeRestraints			*/
/*									*/
/************************************************************************/

ComputeRestraints::ComputeRestraints(ComputeID c, PatchID pid)
  : ComputePatch(c,pid)
{
	reduction = ReductionMgr::Object();
	reduction->Register(REDUCTION_BC_ENERGY);

	SimParameters *simParams = Node::Object()->simParameters;

	//  Get parameters from the SimParameters object
	consExp = simParams->constraintExp;
	
	//****** BEGIN moving constraints changes 
	consMoveOn = simParams->movingConstraintsOn;
	if (consMoveOn) {
	  moveVel = simParams->movingConsVel;
	}
	//****** END moving constraints changes 
	//****** BEGIN rotating constraints changes 
	consRotOn = simParams->rotConstraintsOn;
	if (consRotOn) {
	  rotVel = simParams->rotConsVel;
	  rotAxis = simParams->rotConsAxis;
	  rotPivot = simParams->rotConsPivot;
	}
	//****** END rotating constraints changes 

}
/*			END OF FUNCTION ComputeRestraints		*/

/************************************************************************/
/*									*/
/*			FUNCTION ~ComputeRestraints			*/
/*									*/
/*	This is the destructor for the ComputeRestraints force object.	*/
/*									*/
/************************************************************************/

ComputeRestraints::~ComputeRestraints()

{
	reduction->unRegister(REDUCTION_BC_ENERGY);
}
/*			END OF FUNCTION ~ComputeRestraints		*/

/************************************************************************/
/*									*/
/*				FUNCTION force				*/
/*									*/
/************************************************************************/

void ComputeRestraints::doForce(Position* p, Results* res, AtomProperties* a)

{
	Molecule *molecule = Node::Object()->molecule;
	Real k;			//  Force constant
	Vector refPos;		//  Reference position
	BigReal r, r2; 	//  r=distance between atom position and the
			//  reference position, r2 = r^2
	Vector Rij;	//  vector between current position and reference pos
	BigReal value;	//  Current calculated value

	// aliases to work with old code
	Force *f = res->f[Results::normal];
	BigReal energy = 0;
	BigReal m[9];

	// BEGIN moving and rotating constraint changes ******

	// This version only allows one atom to be moved 
	// and only ALL ref positions to be rotated

	int currentTime = patch->flags.seq;
	if (consRotOn) {
	  vec_rotation_matrix(rotVel * currentTime, rotAxis, m);
	}

	// END moving and rotating constraint changes ******

	  
	for (int localID=0; localID<numAtoms; ++localID)
	{
	  if (molecule->is_atom_constrained(a[localID].id))
	  {
	    molecule->get_cons_params(k, refPos, a[localID].id);

	    // BEGIN moving and rotating constraint changes ******
	    
	    if (consMoveOn) {
	      Rij = refPos + currentTime * moveVel  - p[localID];
	    }
	    else if(consRotOn) {
	      Rij = mat_multiply_vec(refPos - rotPivot, m) + rotPivot 
		- p[localID];
	    }
	    else { // the default case
	     Rij = refPos - p[localID];
	    } 

	    // END moving and rotating constraint changes *******
	    
	    //  Calculate the distance and the distance squared
	    r2 = Rij.length2();
	    r = sqrt(r2);

	    //  Only calculate the energy and the force if the distance is
	    //  non-zero.   Otherwise, you could end up dividing by 0, which
	    //  is bad
	    if (r>0.0)
	    {
	      value=k;
      
	      //  Loop through and multiple k by r consExp times.
	      //  i.e.  calculate kr^e
	      //  I know, I could use pow(), but I don't trust it.
	      for (int k=0; k<consExp; ++k)
	      {
		value *= r;
	      }

	      //  Add to the energy total
	      energy += value;
      
	      //  Now calculate the force, which is ekr^(e-1).  Also, divide
	      //  by another factor of r to normalize the vector before we
	      //  multiple it by the magnitude
	      value *= consExp;
	      value /= r2;
      
	      Rij.mult(value);
      
	      f[localID] += Rij;
	    }
	  }
	}

	reduction->submit(patch->flags.seq, REDUCTION_BC_ENERGY, energy);

}
/*			END OF FUNCTION force				*/

/***************************************************************************
 * RCS INFORMATION:
 *
 *	$RCSfile: ComputeRestraints.C,v $
 *	$Author: sergei $	$Locker:  $		$State: Exp $
 *	$Revision: 1.4 $	$Date: 1998/10/01 00:31:31 $
 *
 ***************************************************************************
 * REVISION HISTORY:
 *
 * $Log: ComputeRestraints.C,v $
 * Revision 1.4  1998/10/01 00:31:31  sergei
 * added rotating restraints feature;
 * changed the moving restraints from only moving one atom to moving all
 * atoms that are restrained. One-atom pulling is available in SMD feature.
 *
 * Revision 1.3  1998/01/05 20:19:04  sergei
 * removed cycling over an array of constrained atoms for moving
 * restraints. Now it's just clean one-atom case.
 *
 * Revision 1.2  1997/08/18 20:16:04  sergei
 * added moving restraint capability with input from config file
 * one atom only
 *
 * Revision 1.1  1997/04/22 04:25:59  jim
 * Added atomic restraints (harmonic constraints) via ComputeRestraints class.
 *
 *
 ***************************************************************************/

