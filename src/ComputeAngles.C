/**
***  Copyright (c) 1995, 1996, 1997, 1998, 1999, 2000 by
***  The Board of Trustees of the University of Illinois.
***  All rights reserved.
**/

/*
   Methods for ComputeAngles.  Main code is for
   loading in the AngleElem information and
   for computing forces and energies for all angles on node's.
   HomePatch(es)
*/

#include "ComputeAngles.h"
#include "Molecule.h"
#include "Parameters.h"
#include "Node.h"
#include "ReductionMgr.h"
#include "Lattice.h"

#include "Debug.h"

#if 0
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
#endif

void AngleElem::getMoleculePointers
    (Molecule* mol, int* count, int32*** byatom, Angle** structarray)
{
  *count = mol->numAngles;
  *byatom = mol->anglesByAtom;
  *structarray = mol->angles;
}

void AngleElem::getParameterPointers(Parameters *p, const AngleValue **v) {
  *v = p->angle_array;
}

void AngleElem::computeForce(BigReal *reduction)
{
  DebugM(3, "::computeForce() localIndex = " << localIndex[0] << " "
               << localIndex[1] << " " << localIndex[2] << endl);

  const Position & pos1 = p[0]->x[localIndex[0]].position;
  const Lattice & lattice = p[0]->p->lattice;
  const Position & pos2 = p[1]->x[localIndex[1]].position;
  const Vector r12 = lattice.delta(pos1,pos2);
  BigReal d12 = r12.length();
  const Position & pos3 = p[2]->x[localIndex[2]].position;
  const Vector r32 = lattice.delta(pos3,pos2);
  BigReal d32 = r32.length();

  BigReal cos_theta = (r12*r32)/(d12*d32);
  //  This code is useless because below we divide by sin_theta!  -JCP
  //  Make sure that the cosine value is acceptable.  With roundoff, you
  //  can get values like 1.0+2e-16, which makes acos puke.  So instead,
  //  just set these kinds of values to exactly 1.0
  // if (cos_theta > 1.0) cos_theta = 1.0;
  // else if (cos_theta < -1.0) cos_theta = -1.0;

  BigReal k = value->k * scale;
  BigReal theta0 = value->theta0;

  //  Get theta
  BigReal theta = acos(cos_theta);

  //  Compare it to the rest angle
  BigReal diff = theta - theta0;

  //  Add the energy from this angle to the total energy
  BigReal energy = k *diff*diff;

  //  Normalize vector r12 and r32
  BigReal d12inv = 1. / d12;
  BigReal d32inv = 1. / d32;

  //  Calculate constant factor 2k(theta-theta0)/sin(theta)
  BigReal sin_theta = sqrt(1.0 - cos_theta*cos_theta);
  diff *= (-2.0* k) / sin_theta;
  BigReal c1 = diff * d12inv;
  BigReal c2 = diff * d32inv;

  //  Calculate the actual forces
  Force force1 = c1*(r12*(d12inv*cos_theta) - r32*d32inv);
  Force force2 = force1;
  Force force3 = c2*(r32*(d32inv*cos_theta) - r12*d12inv);
  force2 += force3;  force2 *= -1;

  //  Check to see if we need to do the Urey-Bradley term
  if (value->k_ub)
  {
	//  Non-zero k_ub value, so calculate the harmonic
	//  potential between the 1-3 atoms
	BigReal k_ub = value->k_ub;
	BigReal r_ub = value->r_ub;
	Vector r13 = r12 - r32;
	BigReal d13 = r13.length();
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
  reduction[virialIndex_XX] += ( force1.x * r12.x + force3.x * r32.x );
  reduction[virialIndex_XY] += ( force1.x * r12.y + force3.x * r32.y );
  reduction[virialIndex_XZ] += ( force1.x * r12.z + force3.x * r32.z );
  reduction[virialIndex_YX] += ( force1.y * r12.x + force3.y * r32.x );
  reduction[virialIndex_YY] += ( force1.y * r12.y + force3.y * r32.y );
  reduction[virialIndex_YZ] += ( force1.y * r12.z + force3.y * r32.z );
  reduction[virialIndex_ZX] += ( force1.z * r12.x + force3.z * r32.x );
  reduction[virialIndex_ZY] += ( force1.z * r12.y + force3.z * r32.y );
  reduction[virialIndex_ZZ] += ( force1.z * r12.z + force3.z * r32.z );
}


void AngleElem::submitReductionData(BigReal *data, SubmitReduction *reduction)
{
  reduction->item(REDUCTION_ANGLE_ENERGY) += data[angleEnergyIndex];
  ADD_TENSOR(reduction,REDUCTION_VIRIAL_NORMAL,data,virialIndex);
}

