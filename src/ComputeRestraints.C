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
	consMoveOn = simParams->movingConstraintsOn;
	if (consMoveOn) {
	  moveVel = simParams->movingConsVel;
	  moveAtom = simParams->movingConsAtom;
	}

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

	// BEGIN moving constraint changes ****************************
	
	// !!!!!!!!! most of the following is a hack, which requires
	// !!!!!!!!! recompilation if the velocity and/or atom numbers
	// !!!!!!!!! change. The slight improvement is with the
	// !!!!!!!!! introduction of moving restraints in the config file,
	// !!!!!!!!! but for now it is only one atom that can be 
	// !!!!!!!!! specified there. In general, it should be a list of
	// !!!!!!!!! atoms.
	// !!!!!!!!! I save the hack just for the heck of it.
	// !!!!!!!!!                8/14/97              Sergei Izrailev

	// The following specifies the velocity of the change of
	// the constraint reference position. For now, it is the
	// same velocity for all reference positions which we want to move.
	// The indices of the atoms for which the reference positions 
	// should be specified here too. 

	int j; // loop counter
	int flag; // temporary variable
	int timestep;

	timestep = patch->flags.seq;

	// put the necessary numbers for the velocity (\AA / timestep) in here

	// The numbers  are for pulling
	// lipid 69 N50 in the pla2 tight structure towards LAM61N50 after
	// 50 ps equilibration in the T-bath + 50 ps free dynamics. 
	// The total displacement should be 
	// 4 \AA in 100 ps.
	//	Vector vel(-1.002556e-05, -3.636444e-06, -3.855208e-05);

	// pulling of  lipid 51 in the loosely coupled pla2-dlpe
	//	Vector vel(1.1938e-05, -5.328e-06, 3.498e-06);
	
	// list the (global) indices of the constrained atoms for which the 
	// reference positions should be moved here, as well as number of 
	// these atoms. NOTE: if the numbers are inconsistent, namd will
	// crash. Atom N in the PDB file corresponds to index N-1 in namd.
	int numToMove = 0; // 0 for equilibration
	
	// the change to take care of input from config file
	if (consMoveOn) { // add 1 to include the atom specified in config file
	  numToMove++;
	}

	int *toMove = new int[numToMove+1]; 
	toMove[0] = moveAtom;  // atom index from the config file
	//   toMove[0] = 3244;  // pla2 loose 
	//   toMove[0] = 4042;  // tight 69_61
	//   toMove[1] = 34;
	//   toMove[2] = 56;
	//   toMove[3] = 67;	

	// END moving constraint changes ****************************


	for (int localID=0; localID<numAtoms; ++localID)
	{
	  if (molecule->is_atom_constrained(a[localID].id))
	  {
	    molecule->get_cons_params(k, refPos, a[localID].id);

	    // BEGIN moving constraint changes ****************************
	    
	    flag = 0;
	    for ( j=0; j<numToMove; j++) {
	      if (a[localID].id == toMove[j]) {
		flag =1;
		break;
	      }
	    }

	    if (flag) { // this reference position is to be moved
	      Rij = refPos + timestep*moveVel - p[localID];
	    }
	    else { // the default case
	      Rij = refPos - p[localID];
	    }

	    // END moving constraint changes ****************************
	    
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
 *	$Revision: 1.2 $	$Date: 1997/08/18 20:16:04 $
 *
 ***************************************************************************
 * REVISION HISTORY:
 *
 * $Log: ComputeRestraints.C,v $
 * Revision 1.2  1997/08/18 20:16:04  sergei
 * added moving restraint capability with input from config file
 * one atom only
 *
 * Revision 1.1  1997/04/22 04:25:59  jim
 * Added atomic restraints (harmonic constraints) via ComputeRestraints class.
 *
 *
 ***************************************************************************/
