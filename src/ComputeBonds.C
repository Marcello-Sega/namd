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

#include "ComputeBonds.h"
#include "Molecule.h"
#include "Parameters.h"
#include "Node.h"
#include "ReductionMgr.h"
#include "Lattice.h"

#include "Debug.h"

void BondElem::loadTuplesForAtom
  (void *voidlist, AtomID atomID, Molecule *molecule)
{
      DebugM(1, "::loadTuplesForAtom - atomID " << atomID << endl );
      UniqueSet<BondElem> &bondList =
                  *( (UniqueSet<BondElem>*) voidlist );

      DebugM(1, "::loadTuplesForAtom - current list size " << bondList.size() << endl );

      /* get list of all bonds for the atom */
      LintList *bonds = molecule->get_bonds_for_atom(atomID);
      DebugM(1, "::loadTuplesForAtom - atomID " << atomID << endl );
      DebugM(1, "::loadTuplesForAtom - bonds->head()" << bonds->head() << endl );

      /* cycle through each bond */
      int bondNum = bonds->head();
      while(bondNum != LIST_EMPTY)
      {
        /* store bond in the list */
        DebugM(1, "::loadTuplesForAtom - loading bond " << bondNum << endl );
        bondList.add(BondElem(molecule->get_bond(bondNum)));
        bondNum = bonds->next();
      }
}

void BondElem::computeForce(BigReal *reduction)
{
  DebugM(1, "::computeForce() localIndex = " << localIndex[0] << " "
               << localIndex[1] << endl);

  Force r12;	// vector between atoms 1,2
  BigReal r;		// Distance between atoms
  BigReal diff;		// difference between theta and theta0
  BigReal energy;	// energy from the bond

  DebugM(3, "::computeForce() -- starting with bond type " << bondType << endl);

  // get the bond information
  Real k, x0;
  Node::Object()->parameters->get_bond_params(&k,&x0,bondType);

  // compute vectors between atoms and their distances
  const Lattice & lattice = p[0]->p->lattice;
  r12 = lattice.delta(p[0]->x[localIndex[0]],p[1]->x[localIndex[1]]);
  r = r12.length();

  //  Compare it to the rest bond
  diff = r - x0;

  //  Add the energy from this bond to the total energy
  energy = k*diff*diff;

  //  Determine the magnitude of the force
  diff *= -2.0*k;

  //  Divide by r to normalize the vector
  diff /= r;

  //  Scale the force vector accordingly
  r12 *= diff;

  //  Now add the forces to each force vector
  p[0]->f[localIndex[0]] += r12;
  p[1]->f[localIndex[1]] -= r12;

  DebugM(3, "::computeForce() -- ending with delta energy " << energy << endl);
  if ( p[0]->patchType == HOME ) reduction[bondEnergyIndex] += energy;
}


void BondElem::registerReductionData(ReductionMgr *reduction)
{
  reduction->Register(REDUCTION_BOND_ENERGY);
}

void BondElem::submitReductionData(BigReal *data, ReductionMgr *reduction, int seq)
{
  reduction->submit(seq, REDUCTION_BOND_ENERGY, data[bondEnergyIndex]);
}

void BondElem::unregisterReductionData(ReductionMgr *reduction)
{
  reduction->unRegister(REDUCTION_BOND_ENERGY);
}

