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

#include "ComputeImpropers.h"
#include "Molecule.h"
#include "Parameters.h"
#include "Node.h"
#include "ReductionMgr.h"
#include "Lattice.h"

#include "Debug.h"

void ImproperElem::loadTuplesForAtom
  (void *voidlist, AtomID atomID, Molecule *molecule)
{
      DebugM(1, "::loadTuplesForAtom - atomID " << atomID << endl );
      UniqueSet<ImproperElem> &improperList =
                  *( (UniqueSet<ImproperElem>*) voidlist );

      DebugM(1, "::loadTuplesForAtom - current list size " << improperList.size() << endl );

      /* get list of all impropers for the atom */
      int *impropers = molecule->get_impropers_for_atom(atomID);
      DebugM(1, "::loadTuplesForAtom - atomID " << atomID << endl );

      /* cycle through each improper */
      int improperNum = *impropers;
      while(improperNum != -1)
      {
        /* store improper in the list */
        DebugM(1, "::loadTuplesForAtom - loading improper " << improperNum << endl );
        improperList.add(ImproperElem(molecule->get_improper(improperNum)));
        improperNum = *(++impropers);
      }
}

void ImproperElem::getMoleculePointers
    (Molecule* mol, int* count, int*** byatom, Improper** structarray)
{
  *count = mol->numImpropers;
  *byatom = mol->impropersByAtom;
  *structarray = mol->impropers;
}

void ImproperElem::computeForce(BigReal *reduction)
{
  DebugM(3, "::computeForce() localIndex = " << localIndex[0] << " "
               << localIndex[1] << " " << localIndex[2] << " " <<
	       localIndex[3] << endl);

  // Vector r12, r23, r34;	// vector between atoms
  Vector A,B,C;		// cross products
  BigReal rA, rB, rC;	// length of vectors A, B, and C
  BigReal energy=0;	// energy from the angle
  BigReal phi;		// angle between the plans
  double cos_phi;	// cos(phi)
  double sin_phi;	// sin(phi)
  Vector dcosdA;	// Derivative d(cos(phi))/dA
  Vector dcosdB;	// Derivative d(cos(phi))/dB
  Vector dsindC;	// Derivative d(sin(phi))/dC
  Vector dsindB;	// Derivative d(sin(phi))/dB
  BigReal K,K1;		// energy constants
  BigReal diff;		// for periodicity
  Force f1,f2,f3;	// force components
  Real k, delta;	// angle information
  int n;		// angle information

  DebugM(3, "::computeForce() -- starting with improper type " << improperType << endl);

  // get the improper information
  int multiplicity = Node::Object()->parameters->get_improper_multiplicity(improperType);

  //  Calculate the vectors between atoms
  const Position & pos0 = p[0]->x[localIndex[0]];
  const Position & pos1 = p[1]->x[localIndex[1]];
  const Position & pos2 = p[2]->x[localIndex[2]];
  const Position & pos3 = p[3]->x[localIndex[3]];
  const Lattice & lattice = p[0]->p->lattice;
  const Vector r12 = lattice.delta(pos0,pos1);
  const Vector r23 = lattice.delta(pos1,pos2);
  const Vector r34 = lattice.delta(pos2,pos3);

  //  Calculate the cross products
  A = cross(r12,r23);
  B = cross(r23,r34);
  C = cross(r23,A);

  //  Calculate the distances
  rA = A.length();
  rB = B.length();
  rC = C.length();

  //  Calculate the sin and cos
  cos_phi = A*B/(rA*rB);
  sin_phi = C*B/(rC*rB);

  //  Normalize B
  rB = 1.0/rB;
  B *= rB;

  errno=0;
  phi= -atan2(sin_phi,cos_phi);
  CHECK_DOMAIN_ATAN(sin_phi,cos_phi);

  if (fabs(sin_phi) > 0.1)
  {
    //  Normalize A
    rA = 1.0/rA;
    A *= rA;
    dcosdA = rA*(cos_phi*A-B);
    dcosdB = rB*(cos_phi*B-A);
  }
  else
  {
    //  Normalize C
    rC = 1.0/rC;
    C *= rC;
    dsindC = rC*(sin_phi*C-B);
    dsindB = rB*(sin_phi*B-C);
  }

  //  Loop through the multiple parameter sets for this
  //  bond.  We will only loop more than once if this
  //  has multiple parameter sets from Charmm22
  for (int mult_num=0; mult_num<multiplicity; mult_num++)
  {
    /* get angle information */
    Node::Object()->parameters->get_improper_params(&k,&n,&delta,improperType,mult_num);

    //  Calculate the energy
    if (n)
    {
      //  Periodicity is greater than 0, so use cos form
      K = k*(1+cos(n*phi + delta));
      K1 = -n*k*sin(n*phi + delta);
    }
    else
    {
      //  Periodicity is 0, so just use the harmonic form
      diff = phi-delta;
      if (diff < -PI)           diff += TWOPI;
      else if (diff > PI)       diff -= TWOPI;

      K = k*diff*diff;
      K1 = 2.0*k*diff;
    }

    //  Add the energy from this improper to the total energy
    energy += K;

    //  Next, we want to calculate the forces.  In order
    //  to do that, we first need to figure out whether the
    //  sin or cos form will be more stable.  For this,
    //  just look at the value of phi
    if (fabs(sin_phi) > 0.1)
    {
      //  use the sin version to avoid 1/cos terms
      K1 = K1/sin_phi;

      f1.x += K1*(r23.y*dcosdA.z - r23.z*dcosdA.y);
      f1.y += K1*(r23.z*dcosdA.x - r23.x*dcosdA.z);
      f1.z += K1*(r23.x*dcosdA.y - r23.y*dcosdA.x);

      f3.x += K1*(r23.z*dcosdB.y - r23.y*dcosdB.z);
      f3.y += K1*(r23.x*dcosdB.z - r23.z*dcosdB.x);
      f3.z += K1*(r23.y*dcosdB.x - r23.x*dcosdB.y);

      f2.x += K1*(r12.z*dcosdA.y - r12.y*dcosdA.z
               + r34.y*dcosdB.z - r34.z*dcosdB.y);
      f2.y += K1*(r12.x*dcosdA.z - r12.z*dcosdA.x
               + r34.z*dcosdB.x - r34.x*dcosdB.z);
      f2.z += K1*(r12.y*dcosdA.x - r12.x*dcosdA.y
             + r34.x*dcosdB.y - r34.y*dcosdB.x);
    }
    else
    {
      //  This angle is closer to 0 or 180 than it is to
      //  90, so use the cos version to avoid 1/sin terms
      K1 = -K1/cos_phi;

      f1.x += K1*((r23.y*r23.y + r23.z*r23.z)*dsindC.x
                - r23.x*r23.y*dsindC.y
                - r23.x*r23.z*dsindC.z);
      f1.y += K1*((r23.z*r23.z + r23.x*r23.x)*dsindC.y
                - r23.y*r23.z*dsindC.z
                - r23.y*r23.x*dsindC.x);
      f1.z += K1*((r23.x*r23.x + r23.y*r23.y)*dsindC.z
                - r23.z*r23.x*dsindC.x
                - r23.z*r23.y*dsindC.y);

      f3 += cross(K1,dsindB,r23);

      f2.x += K1*(-(r23.y*r12.y + r23.z*r12.z)*dsindC.x
             +(2.0*r23.x*r12.y - r12.x*r23.y)*dsindC.y
             +(2.0*r23.x*r12.z - r12.x*r23.z)*dsindC.z
             +dsindB.z*r34.y - dsindB.y*r34.z);
      f2.y += K1*(-(r23.z*r12.z + r23.x*r12.x)*dsindC.y
             +(2.0*r23.y*r12.z - r12.y*r23.z)*dsindC.z
             +(2.0*r23.y*r12.x - r12.y*r23.x)*dsindC.x
             +dsindB.x*r34.z - dsindB.z*r34.x);
      f2.z += K1*(-(r23.x*r12.x + r23.y*r12.y)*dsindC.z
             +(2.0*r23.z*r12.x - r12.z*r23.x)*dsindC.x
             +(2.0*r23.z*r12.y - r12.z*r23.y)*dsindC.y
             +dsindB.y*r34.x - dsindB.x*r34.y);
    }
  } /* for multiplicity */

  /* store the forces */
  p[0]->f[localIndex[0]] += f1;
  p[1]->f[localIndex[1]] += f2 - f1;
  p[2]->f[localIndex[2]] += f3 - f2;
  p[3]->f[localIndex[3]] += -f3;

  DebugM(3, "::computeForce() -- ending with delta energy " << energy << endl);
  reduction[improperEnergyIndex] += energy;
  reduction[virialXIndex] += ( f1.x * r12.x + f2.x * r23.x + f3.x * r34.x );
  reduction[virialYIndex] += ( f1.y * r12.y + f2.y * r23.y + f3.y * r34.y );
  reduction[virialZIndex] += ( f1.z * r12.z + f2.z * r23.z + f3.z * r34.z );
}


void ImproperElem::registerReductionData(ReductionMgr *reduction)
{
  reduction->Register(REDUCTION_IMPROPER_ENERGY);
  reduction->Register(REDUCTION_VIRIAL_NORMAL_X);
  reduction->Register(REDUCTION_VIRIAL_NORMAL_Y);
  reduction->Register(REDUCTION_VIRIAL_NORMAL_Z);
}

void ImproperElem::submitReductionData(BigReal *data, ReductionMgr *reduction, int seq)
{
  reduction->submit(seq, REDUCTION_IMPROPER_ENERGY, data[improperEnergyIndex]);
  reduction->submit(seq, REDUCTION_VIRIAL_NORMAL_X, data[virialXIndex]);
  reduction->submit(seq, REDUCTION_VIRIAL_NORMAL_Y, data[virialYIndex]);
  reduction->submit(seq, REDUCTION_VIRIAL_NORMAL_Z, data[virialZIndex]);
}

void ImproperElem::unregisterReductionData(ReductionMgr *reduction)
{
  reduction->unRegister(REDUCTION_IMPROPER_ENERGY);
  reduction->unRegister(REDUCTION_VIRIAL_NORMAL_X);
  reduction->unRegister(REDUCTION_VIRIAL_NORMAL_Y);
  reduction->unRegister(REDUCTION_VIRIAL_NORMAL_Z);
}



/***************************************************************************
 * RCS INFORMATION:
 *
 *	$RCSfile $
 *	$Author $	$Locker:  $		$State: Exp $
 *	$Revision: 1.1010 $	$Date: 1999/01/06 00:56:22 $
 *
 ***************************************************************************
 * REVISION HISTORY:
 *
 * $Log: ComputeImpropers.C,v $
 * Revision 1.1010  1999/01/06 00:56:22  jim
 * All compute objects except DPMTA now return diagonal of virial tensor.
 *
 * Revision 1.1009  1998/08/03 22:09:33  brunner
 * Cleared errno before trig functions
 *
 * Revision 1.1008  1998/06/18 14:48:02  jim
 * Split virial into NORMAL, NBOND, and SLOW parts to match force classes.
 *
 * Revision 1.1007  1997/10/17 17:16:47  jim
 * Switched from hash tables to checklists, eliminated special exclusion code.
 *
 * Revision 1.1006  1997/09/28 22:36:50  jim
 * Modified tuple-based computations to not duplicate calculations and
 * only require "upstream" proxies.
 *
 * Revision 1.1005  97/08/18  05:02:54  05:02:54  jim (Jim Phillips)
 * Fixed bugs related to multiple dihedrals and exclude 1-2, 1-3, or 1-4.
 * 
 * Revision 1.1004  1997/03/19 11:54:09  ari
 * Add Broadcast mechanism.
 * Fixed RCS Log entries on files that did not have Log entries.
 * Added some register variables to Molecule and ComputeNonbondedExcl.C
 *
 ***************************************************************************/
