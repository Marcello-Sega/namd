/***************************************************************************/
/*          (C) Copyright 1996,1997 The Board of Trustees of th            */
/*                          University of Illinois                         */
/*                           All Rights Reserved                           */
/***************************************************************************/
/***************************************************************************
 * DESCRIPTION: Methods for ComputeAngles.  Main code is for
 *		loading in the AngleElem information and
 *		for computing forces and energies for all angles on node's.
 *		HomePatch(es)
 *
 ***************************************************************************/

#include "ComputeAngles.h"
#include "Molecule.h"
#include "Parameters.h"
#include "Node.h"
#include "ReductionMgr.h"
#include "Lattice.h"

#include "Debug.h"

void AngleElem::loadTuplesForAtom
  (void *voidlist, AtomID atomID, Molecule *molecule)
{
      DebugM(1, "::loadTuplesForAtom - atomID " << atomID << endl );
      UniqueSet<AngleElem> &angleList =
                  *( (UniqueSet<AngleElem>*) voidlist );

      DebugM(1, "::loadTuplesForAtom - current list size " << angleList.size() << endl );

      /* get list of all angles for the atom */
      int *angles = molecule->get_angles_for_atom(atomID);
      DebugM(1, "::loadTuplesForAtom - atomID " << atomID << endl );

      /* cycle through each angle */
      int angleNum = *angles;
      while(angleNum != -1)
      {
        /* store angle in the list */
        DebugM(1, "::loadTuplesForAtom - loading angle " << angleNum << endl );
        angleList.add(AngleElem(molecule->get_angle(angleNum)));
        angleNum = *(++angles);
      }
}

void AngleElem::getMoleculePointers
    (Molecule* mol, int* count, int*** byatom, Angle** structarray)
{
  *count = mol->numAngles;
  *byatom = mol->anglesByAtom;
  *structarray = mol->angles;
}

void AngleElem::computeForce(BigReal *reduction)
{
  DebugM(3, "::computeForce() localIndex = " << localIndex[0] << " "
               << localIndex[1] << " " << localIndex[2] << endl);

  const Position & pos1 = p[0]->x[localIndex[0]];
  const Position & pos2 = p[1]->x[localIndex[1]];
  const Position & pos3 = p[2]->x[localIndex[2]];
  Force force1;
  Force force2;
  Force force3;
  const Lattice & lattice = p[0]->p->lattice;

  BigReal d12, d32, d13;	// distances between atoms
  BigReal theta;	// theta
  BigReal cos_theta;	// cos(theta)
  BigReal sin_theta;	// sin(theta)
  register BigReal diff;		// difference between theta and theta0
  BigReal c1,c2;	// constant factors involved in force
  BigReal energy;	// energy from the angle

  DebugM(3, "::computeForce() -- starting with angle type " << angleType << endl);

  // get the angle information
  Real k, theta0, k_ub, r_ub;
  Node::Object()->parameters->get_angle_params(&k,&theta0,&k_ub,&r_ub,angleType);

  // compute vectors between atoms and their distances
  const Vector r12 = lattice.delta(pos1,pos2);
  const Vector r32 = lattice.delta(pos3,pos2);

  d12 = r12.length();
  d32 = r32.length();

  DebugM(3, "::computeForce() d12 = " << d12 << " d32 = " << d32 << endl);

  //  Make sure that the cosine value is acceptable.  With roundoff, you
  //  can get values like 1.0+2e-16, which makes acos puke.  So instead,
  //  just set these kinds of values to exactly 1.0
  cos_theta = (r12*r32)/(d12*d32);
  if (cos_theta > 1.0) cos_theta = 1.0;
  else if (cos_theta < -1.0) cos_theta = -1.0;
  sin_theta = sqrt(1.0 - cos_theta*cos_theta);

  //  Get theta
  errno=0;
  theta = acos(cos_theta);
  CHECK_DOMAIN_ACOS(cos_theta);

  //  Compare it to the rest angle
  diff = theta - theta0;

  DebugM(3, "::computeForce() theta = " << theta << " theta0 = " << theta0 << endl);

  //  Add the energy from this angle to the total energy
  energy = k *diff*diff;

  //  Normalize vector r12 and r32
  BigReal d12inv = 1. / d12;
  BigReal d32inv = 1. / d32;

  //  Calculate constant factor 2k(theta-theta0)/sin(theta)
  diff *= (-2.0* k) / sin_theta;
  c1 = diff * d12inv;
  c2 = diff * d32inv;

  //  Calculate the actual forces
  Force f12 = c1*(r12*(d12inv*cos_theta) - r32*d32inv);
  force1 += f12; force2 -= f12;
  Force f32 = c2*(r32*(d32inv*cos_theta) - r12*d12inv);
  force3 += f32; force2 -= f32;

  //  Check to see if we need to do the Urey-Bradley term
  //  Forces are used only when atom1 or atom3 are local
  if (k_ub > 0.0)
  {
	//  Non-zero k_ub value, so calculate the harmonic
	//  potential between the 1-3 atoms
	Vector r13 = r12 - r32;
	d13 = r13.length();
	diff = d13- r_ub;

	energy += k_ub *diff*diff;

	diff *= -2.0*k_ub / d13;
	r13 *= diff;

	force1 += r13;
	force3 -= r13;
  }

  p[0]->f[localIndex[0]] += force1;
  p[1]->f[localIndex[1]] += force2;
  p[2]->f[localIndex[2]] += force3;

  DebugM(3, "::computeForce() -- ending with delta energy " << energy << endl);
  reduction[angleEnergyIndex] += energy;
  reduction[virialXIndex] += ( r12.x * force1.x + r32.x * force3.x );
  reduction[virialYIndex] += ( r12.y * force1.y + r32.y * force3.y );
  reduction[virialZIndex] += ( r12.z * force1.z + r32.z * force3.z );
}


void AngleElem::registerReductionData(ReductionMgr *reduction)
{
  reduction->Register(REDUCTION_ANGLE_ENERGY);
  reduction->Register(REDUCTION_VIRIAL_NORMAL_X);
  reduction->Register(REDUCTION_VIRIAL_NORMAL_Y);
  reduction->Register(REDUCTION_VIRIAL_NORMAL_Z);
}

void AngleElem::submitReductionData(BigReal *data, ReductionMgr *reduction, int seq)
{
  reduction->submit(seq, REDUCTION_ANGLE_ENERGY, data[angleEnergyIndex]);
  reduction->submit(seq, REDUCTION_VIRIAL_NORMAL_X, data[virialXIndex]);
  reduction->submit(seq, REDUCTION_VIRIAL_NORMAL_Y, data[virialYIndex]);
  reduction->submit(seq, REDUCTION_VIRIAL_NORMAL_Z, data[virialZIndex]);
}

void AngleElem::unregisterReductionData(ReductionMgr *reduction)
{
  reduction->unRegister(REDUCTION_ANGLE_ENERGY);
  reduction->unRegister(REDUCTION_VIRIAL_NORMAL_X);
  reduction->unRegister(REDUCTION_VIRIAL_NORMAL_Y);
  reduction->unRegister(REDUCTION_VIRIAL_NORMAL_Z);
}

/***************************************************************************
 * RCS INFORMATION:
 *
 *	$RCSfile $
 *	$Author $	$Locker:  $		$State: Exp $
 *	$Revision: 1.1010 $	$Date: 1999/01/06 00:56:19 $
 *
 ***************************************************************************
 * REVISION HISTORY:
 *
 * $Log: ComputeAngles.C,v $
 * Revision 1.1010  1999/01/06 00:56:19  jim
 * All compute objects except DPMTA now return diagonal of virial tensor.
 *
 * Revision 1.1009  1998/08/03 22:09:31  brunner
 * Cleared errno before trig functions
 *
 * Revision 1.1008  1998/06/18 14:47:59  jim
 * Split virial into NORMAL, NBOND, and SLOW parts to match force classes.
 *
 * Revision 1.1007  1997/10/17 17:16:42  jim
 * Switched from hash tables to checklists, eliminated special exclusion code.
 *
 * Revision 1.1006  1997/09/28 22:36:47  jim
 * Modified tuple-based computations to not duplicate calculations and
 * only require "upstream" proxies.
 *
 * Revision 1.1005  1997/03/20 23:53:27  ari
 * Some changes for comments. Copyright date additions.
 * Hooks for base level update of Compute objects from ComputeMap
 * by ComputeMgr.  Useful for new compute migration functionality.
 *
 * Revision 1.1004  1997/03/19 11:54:02  ari
 * Add Broadcast mechanism.
 * Fixed RCS Log entries on files that did not have Log entries.
 * Added some register variables to Molecule and ComputeNonbondedExcl.C
 *
 ***************************************************************************/
