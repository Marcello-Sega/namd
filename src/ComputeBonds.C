/**
***  Copyright (c) 1995, 1996, 1997, 1998, 1999, 2000 by
***  The Board of Trustees of the University of Illinois.
***  All rights reserved.
**/

#include "ComputeBonds.h"
#include "Molecule.h"
#include "Parameters.h"
#include "Node.h"
#include "ReductionMgr.h"
#include "Lattice.h"

#include "Debug.h"

#if 0
void BondElem::loadTuplesForAtom
  (void *voidlist, AtomID atomID, Molecule *molecule)
{
      DebugM(1, "::loadTuplesForAtom - atomID " << atomID << endl );
      UniqueSet<BondElem> &bondList =
                  *( (UniqueSet<BondElem>*) voidlist );

      DebugM(1, "::loadTuplesForAtom - current list size " << bondList.size() << endl );

      /* get list of all bonds for the atom */
      int *bonds = molecule->get_bonds_for_atom(atomID);
      DebugM(1, "::loadTuplesForAtom - atomID " << atomID << endl );

      /* cycle through each bond */
      int bondNum = *bonds;
      while(bondNum != -1)
      {
        /* store bond in the list */
        DebugM(1, "::loadTuplesForAtom - loading bond " << bondNum << endl );
        bondList.add(BondElem(molecule->get_bond(bondNum)));
        bondNum = *(++bonds);
      }
}
#endif

void BondElem::getMoleculePointers
    (Molecule* mol, int* count, int*** byatom, Bond** structarray)
{
  *count = mol->numBonds;
  *byatom = mol->bondsByAtom;
  *structarray = mol->bonds;
}

void BondElem::getParameterPointers(Parameters *p, const BondValue **v) {
  *v = p->bond_array;
}

void BondElem::computeForce(BigReal *reduction)
{
  DebugM(1, "::computeForce() localIndex = " << localIndex[0] << " "
               << localIndex[1] << endl);

  BigReal r;		// Distance between atoms
  BigReal diff;		// difference between theta and theta0
  BigReal energy;	// energy from the bond

  DebugM(3, "::computeForce() -- starting with bond type " << bondType << endl);

  // get the bond information
  Real k = value->k;
  Real x0 = value->x0;

  // compute vectors between atoms and their distances
  const Lattice & lattice = p[0]->p->lattice;
  const Vector r12 = lattice.delta(p[0]->x[localIndex[0]].position,
					p[1]->x[localIndex[1]].position);
  r = r12.length();

  //  Compare it to the rest bond
  diff = r - x0;

  //  Add the energy from this bond to the total energy
  energy = k*diff*diff;

  //  Determine the magnitude of the force
  diff *= -2.0*k;

  //  Scale the force vector accordingly
  const Force f12 = r12 * (diff/r);

  //  Now add the forces to each force vector
  p[0]->f[localIndex[0]] += f12;
  p[1]->f[localIndex[1]] -= f12;

  DebugM(3, "::computeForce() -- ending with delta energy " << energy << endl);
  reduction[bondEnergyIndex] += energy;
  reduction[virialIndex_XX] += f12.x * r12.x;
  reduction[virialIndex_XY] += f12.x * r12.y;
  reduction[virialIndex_XZ] += f12.x * r12.z;
  reduction[virialIndex_YX] += f12.y * r12.x;
  reduction[virialIndex_YY] += f12.y * r12.y;
  reduction[virialIndex_YZ] += f12.y * r12.z;
  reduction[virialIndex_ZX] += f12.z * r12.x;
  reduction[virialIndex_ZY] += f12.z * r12.y;
  reduction[virialIndex_ZZ] += f12.z * r12.z;
}


void BondElem::submitReductionData(BigReal *data, SubmitReduction *reduction)
{
  reduction->item(REDUCTION_BOND_ENERGY) += data[bondEnergyIndex];
  ADD_TENSOR(reduction,REDUCTION_VIRIAL_NORMAL,data,virialIndex);
}

