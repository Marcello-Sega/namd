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

#include "ComputeDihedrals.h"
#include "Molecule.h"
#include "Parameters.h"
#include "Node.h"
#include "Debug.h"
#include "ReductionMgr.h"
#include "Lattice.h"

void DihedralElem::loadTuplesForAtom
  (void *voidlist, AtomID atomID, Molecule *molecule)
{
      DebugM(1, "::loadTuplesForAtom - atomID " << atomID << endl );
      UniqueSet<DihedralElem> &dihedralList =
                  *( (UniqueSet<DihedralElem>*) voidlist );

      DebugM(1, "::loadTuplesForAtom - current list size " << dihedralList.size() << endl );

      /* get list of all dihedrals for the atom */
      LintList *dihedrals = molecule->get_dihedrals_for_atom(atomID);
      DebugM(1, "::loadTuplesForAtom - atomID " << atomID << endl );
      DebugM(1, "::loadTuplesForAtom - dihedrals->head()" << dihedrals->head() << endl );

      /* cycle through each dihedral */
      int dihedralNum = dihedrals->head();
      while(dihedralNum != LIST_EMPTY)
      {
        /* store dihedral in the list */
        DebugM(1, "::loadTuplesForAtom - loading dihedral " << dihedralNum << endl );
        dihedralList.add(DihedralElem(molecule->get_dihedral(dihedralNum)));
        dihedralNum = dihedrals->next();
      }
}

void DihedralElem::computeForce(BigReal *reduction)
{
  DebugM(3, "::computeForce() localIndex = " << localIndex[0] << " "
               << localIndex[1] << " " << localIndex[2] << endl);

  Vector r12, r23, r34;	// vector between atoms
  Vector A,B,C;		// cross products
  BigReal rA, rB, rC;	// length of vectors A, B, and C
  BigReal energy=0;	// energy from the angle
  BigReal phi;		// angle between the plans
  BigReal cos_phi;	// cos(phi)
  BigReal sin_phi;	// sin(phi)
  Vector dcosdA;	// Derivative d(cos(phi))/dA
  Vector dcosdB;	// Derivative d(cos(phi))/dB
  Vector dsindC;	// Derivative d(sin(phi))/dC
  Vector dsindB;	// Derivative d(sin(phi))/dB
  BigReal K,K1;		// energy constants
  BigReal diff;		// for periodicity
  Force f1,f2,f3;	// force components
  Real k, delta;	// angle information
  int n;		// angle information

  DebugM(3, "::computeForce() -- starting with dihedral type " << dihedralType << endl);

  // get the dihedral information
  int multiplicity = Node::Object()->parameters->get_dihedral_multiplicity(dihedralType);

  //  Calculate the vectors between atoms
  const Position & pos0 = p[0]->x[localIndex[0]];
  const Position & pos1 = p[1]->x[localIndex[1]];
  const Position & pos2 = p[2]->x[localIndex[2]];
  const Position & pos3 = p[3]->x[localIndex[3]];
  const Lattice & lattice = p[0]->p->lattice;
  r12 = lattice.delta(pos0,pos1);
  r23 = lattice.delta(pos1,pos2);
  r34 = lattice.delta(pos2,pos3);

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
    Node::Object()->parameters->get_dihedral_params(&k,&n,&delta,dihedralType,mult_num);

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

    //  Add the energy from this dihedral to the total energy
    energy += K;

    //  Next, we want to calculate the forces.  In order
    //  to do that, we first need to figure out whether the
    //  sin or cos form will be more stable.  For this,
    //  just look at the value of phi
    if (fabs(sin_phi) > 0.1)
    {
      //  use the sin version to avoid 1/cos terms
      K1 = K1/sin_phi;

      #if 0
      /** this is not a cross product...
      "Explore" shows it as a very complex
      function.  Double checked with docs and with Robert. -NK **/
      f1 = cross(K1,r23,dcosdA);
      f1.y = K1*(r23.z*dcosdA.x - r23.x*dcosdA.z);


      f3 = cross(K1,dcosdB,r23);
      f3.y = K1*(r23.x*dcosdB.z - r23.z*dcosdB.x);

      f2 = K1*(cross(dcosdA,r12) + cross(r34,dcosdB));
      f2.y = K1*(r12.x*dcosdA.z - r12.z*dcosdA.x
               + r34.z*dcosdB.x - r34.x*dcosdB.z);
      #else
      f1.x = K1*(r23.y*dcosdA.z - r23.z*dcosdA.y);
      f1.y = K1*(r23.z*dcosdA.x - r23.x*dcosdA.z);
      f1.z = K1*(r23.x*dcosdA.y - r23.y*dcosdA.x);

      f3.x = K1*(r23.z*dcosdB.y - r23.y*dcosdB.z);
      f3.y = K1*(r23.x*dcosdB.z - r23.z*dcosdB.x);
      f3.z = K1*(r23.y*dcosdB.x - r23.x*dcosdB.y);

      f2.x = K1*(r12.z*dcosdA.y - r12.y*dcosdA.z
               + r34.y*dcosdB.z - r34.z*dcosdB.y);
      f2.y = K1*(r12.x*dcosdA.z - r12.z*dcosdA.x
               + r34.z*dcosdB.x - r34.x*dcosdB.z);
      f2.z = K1*(r12.y*dcosdA.x - r12.x*dcosdA.y
             + r34.x*dcosdB.y - r34.y*dcosdB.x);
      #endif
    }
    else
    {
      //  This angle is closer to 0 or 180 than it is to
      //  90, so use the cos version to avoid 1/sin terms
      K1 = -K1/cos_phi;

      f1.x = K1*((r23.y*r23.y + r23.z*r23.z)*dsindC.x
                - r23.x*r23.y*dsindC.y
                - r23.x*r23.z*dsindC.z);
      f1.y = K1*((r23.z*r23.z + r23.x*r23.x)*dsindC.y
                - r23.y*r23.z*dsindC.z
                - r23.y*r23.x*dsindC.x);
      f1.z = K1*((r23.x*r23.x + r23.y*r23.y)*dsindC.z
                - r23.z*r23.x*dsindC.x
                - r23.z*r23.y*dsindC.y);

      f3 = cross(K1,dsindB,r23);

      f2.x = K1*(-(r23.y*r12.y + r23.z*r12.z)*dsindC.x
             +(2.0*r23.x*r12.y - r12.x*r23.y)*dsindC.y
             +(2.0*r23.x*r12.z - r12.x*r23.z)*dsindC.z
             +dsindB.z*r34.y - dsindB.y*r34.z);
      f2.y = K1*(-(r23.z*r12.z + r23.x*r12.x)*dsindC.y
             +(2.0*r23.y*r12.z - r12.y*r23.z)*dsindC.z
             +(2.0*r23.y*r12.x - r12.y*r23.x)*dsindC.x
             +dsindB.x*r34.z - dsindB.z*r34.x);
      f2.z = K1*(-(r23.x*r12.x + r23.y*r12.y)*dsindC.z
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
  if ( p[0]->patchType == HOME ) reduction[dihedralEnergyIndex] += energy;
}


void DihedralElem::registerReductionData(ReductionMgr *reduction)
{
  reduction->Register(REDUCTION_DIHEDRAL_ENERGY);
}

void DihedralElem::submitReductionData(BigReal *data, ReductionMgr *reduction, int seq)
{
  reduction->submit(seq, REDUCTION_DIHEDRAL_ENERGY, data[dihedralEnergyIndex]);
}

void DihedralElem::unregisterReductionData(ReductionMgr *reduction)
{
  reduction->unRegister(REDUCTION_DIHEDRAL_ENERGY);
}

