/***************************************************************************/
/*                                                                         */
/*              (C) Copyright 1996 The Board of Trustees of the            */
/*                          University of Illinois                         */
/*                           All Rights Reserved                           */
/*									   */
/***************************************************************************/

/***************************************************************************
 * DESCRIPTION:
 *
 ***************************************************************************/

#include "ComputeAngles.h"
#include "Molecule.h"
#include "Parameters.h"
#include "Node.h"


#define DEBUGM
#include "Debug.h"

void AngleElem::addTuplesForAtom
  (void *voidlist, AtomID atomID, Molecule *molecule)
{
      DebugM(1, "::addTuplesForAtom - atomID " << atomID << endl );
      UniqueSortedArray<AngleElem> &angleList =
                  *( (UniqueSortedArray<AngleElem>*) voidlist );

      DebugM(1, "::addTuplesForAtom - current list size " << angleList.size() << endl );

      /* get list of all angles for the atom */
      LintList *angles = molecule->get_angles_for_atom(atomID);
      DebugM(1, "::addTuplesForAtom - atomID " << atomID << endl );
      DebugM(1, "::addTuplesForAtom - angles->head()" << angles->head() << endl );

      /* cycle through each angle */
      int angleNum = angles->head();
      while(angleNum != LIST_EMPTY)
      {
        /* store angle in the list */
        DebugM(1, "::addTuplesForAtom - adding angle " << angleNum << endl );
        angleList.add(AngleElem(molecule->get_angle(angleNum)));
        angleNum = angles->next();
      }
}

BigReal AngleElem::computeForce(void)
{
  DebugM(3, "::computeForce() localIndex = " << localIndex[0] << " "
               << localIndex[1] << " " << localIndex[2] << endl);

  const Position & pos1 = p[0]->x[localIndex[0]];
  const Position & pos2 = p[1]->x[localIndex[1]];
  const Position & pos3 = p[2]->x[localIndex[2]];
  Force & force1 = p[0]->f[localIndex[0]];
  Force & force2 = p[1]->f[localIndex[1]];
  Force & force3 = p[2]->f[localIndex[2]];

  Vector r12, r32, r13;	// vector between atoms 1,2 and 3,2
  BigReal d12, d32, d13;	// distances between atoms
  BigReal theta;	// theta
  BigReal cos_theta;	// cos(theta)
  BigReal sin_theta;	// sin(theta)
  BigReal diff;		// difference between theta and theta0
  BigReal c1,c2;	// constant factors involved in force
  BigReal energy;	// energy from the angle

  DebugM(3, "::computeForce() -- starting with angle type " << angleType << endl);

  // get the angle information
  Real k, theta0, k_ub, r_ub;
  Node::Object()->parameters->get_angle_params(&k,&theta0,&k_ub,&r_ub,angleType);

  // compute vectors between atoms and their distances
  r12 = pos1-pos2;
  r32 = pos3-pos1;

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
  theta = acos(cos_theta);
  CHECK_DOMAIN_ACOS(cos_theta);

  //  Compare it to the rest angle
  diff = theta - theta0;

  DebugM(3, "::computeForce() theta = " << theta << " theta0 = " << theta0 << endl);

  //  Add the energy from this angle to the total energy
  energy = k *diff*diff;

  //  Normalize vector r12 and r32
  r12 /= d12;
  r32 /= d32;

  //  Calculate constant factor 2k(theta-theta0)/sin(theta)
  diff *= (-2.0* k) / sin_theta;
  c1 = diff/d12;
  c2 = diff/d32;

  //  Calculate the actual forces
  r13 = c1*(r12*cos_theta - r32);
  force1 += r13; force2 -= r13;
  r13 = c2*(r32*cos_theta - r12);
  force3 += r13; force2 -= r13;

  //  Check to see if we need to do the Urey-Bradley term
  //  Forces are used only when atom1 or atom3 are local
  if (k_ub > 0.0)
  {
	//  Non-zero k_ub value, so calculate the harmonic
	//  potential between the 1-3 atoms
	r13 = pos1-pos3;
	d13 = r13.length();
	diff = d13- r_ub;

	energy += k_ub *diff*diff;

	diff *= -2.0*k_ub / d13;
	r13 *= diff;

	force1 += r13;
	force3 -= r13;
  }

  DebugM(3, "::computeForce() -- ending with delta energy " << energy << endl);
  return(energy);
}

