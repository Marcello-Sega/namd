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

	for (int localID=0; localID<numAtoms; ++localID)
	{
	  if (molecule->is_atom_constrained(a[localID].id))
	  {
	    molecule->get_cons_params(k, refPos, a[localID].id);
	    Rij = refPos - p[localID];

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
 *	$Author: jim $	$Locker:  $		$State: Exp $
 *	$Revision: 1.1 $	$Date: 1997/04/22 04:25:59 $
 *
 ***************************************************************************
 * REVISION HISTORY:
 *
 * $Log: ComputeRestraints.C,v $
 * Revision 1.1  1997/04/22 04:25:59  jim
 * Added atomic restraints (harmonic constraints) via ComputeRestraints class.
 *
 *
 ***************************************************************************/
