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

BigReal angleForce (
		const Position pos1, const Position pos2, const Position pos3,
		Force *force1, Force *force2, Force *force3,
		const Index angleType);

BigReal AngleElem::computeForce(void)
{
    return
    angleForce(p[0]->x[localIndex[0]],
	       p[1]->x[localIndex[1]],
	       p[2]->x[localIndex[2]],
	       p[0]->f+localIndex[0],
	       p[1]->f+localIndex[1],
	       p[2]->f+localIndex[2],
	       angleType);
}



/************************************************************************ 
/*									*/
/*			FUNCTION angleForce				*/
/*									*/
/*   INPUTS:								*/
/*      pos1,2,3 - position of atoms					*/
/*      angleType - desired angle (for k, k_ub, r,ub, theta0)		*/
/*									*/
/*   OUTPUTS:								*/
/*      force1,2,3 - forces to be added to the atoms			*/
/*		These forces are not initialized here!  Just added to..	*/
/*	returns - energy from the angle					*/
/*									*/
/************************************************************************/
BigReal angleForce (
		const Position pos1, const Position pos2, const Position pos3,
		Force *force1, Force *force2, Force *force3,
		const Index angleType)
{
  Vector r12, r32, r13;	// vector between atoms 1,2 and 3,2
  BigReal d12, d32, d13;	// distances between atoms
  BigReal theta;	// theta
  BigReal cos_theta;	// cos(theta)
  BigReal sin_theta;	// sin(theta)
  BigReal diff;		// difference between theta and theta0
  BigReal c1,c2;	// constant factors involved in force
  BigReal energy;	// energy from the angle

  CPrintf("ComputeAngles::angleForce() -- starting wiht angle type %d\n",(int)angleType);

  // get the angle information
  Real k, theta0, k_ub, r_ub;
  Node::Object()->parameters->get_angle_params(&k,&theta0,&k_ub,&r_ub,angleType);

  // compute vectors between atoms and their distances
  r12 = pos1-pos2;
  r32 = pos3-pos1;

  d12 = r12.length();
  d32 = r32.length();

  //  Make sure that the cosine value is acceptable.  With roundoff, you
  //  can get values like 1.0+2e-16, which makes acos puke.  So instead,
  //  just set these kinds of values to exactly 1.0
  cos_theta = (r12*r32)/(d12*d32);
  if (cos_theta > 1.0) cos_theta = 1.0;
  else if (cos_theta < -1.0) cos_theta = -1.0;
  sin_theta = sqrt(1.0 - cos_theta*cos_theta);

  //  Get theta
  theta = acos(cos_theta);
  CHECK_DOMAIN();

  //  Compare it to the rest angle
  diff = theta - theta0;

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
  *force1 += c1*(r12*cos_theta - r32);;
  *force3 += c2*(r32*cos_theta - r12);;
  *force2 -= (*force1+*force3);

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

	*force1 += r13;
	*force3 -= r13;
  }

  CPrintf("ComputeAngles::angleForce() -- ending with delta energy %f\n",(float)energy);
  return(energy);
}

