/***************************************************************************/
/*        (C) Copyright 1996,1997 The Board of Trustees of the             */
/*                          University of Illinois                         */
/*                           All Rights Reserved                           */
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

void BondElem::getMoleculePointers
    (Molecule* mol, int* count, int*** byatom, Bond** structarray)
{
  *count = mol->numBonds;
  *byatom = mol->bondsByAtom;
  *structarray = mol->bonds;
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
  Real k, x0;
  Node::Object()->parameters->get_bond_params(&k,&x0,bondType);

  // compute vectors between atoms and their distances
  const Lattice & lattice = p[0]->p->lattice;
  const Vector r12 = lattice.delta(p[0]->x[localIndex[0]],p[1]->x[localIndex[1]]);
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
  reduction[virialXIndex] += r12.x * f12.x;
  reduction[virialYIndex] += r12.y * f12.y;
  reduction[virialZIndex] += r12.z * f12.z;
}


void BondElem::submitReductionData(BigReal *data, SubmitReduction *reduction, int seq)
{
  reduction->item(REDUCTION_BOND_ENERGY) += data[bondEnergyIndex];
  reduction->item(REDUCTION_VIRIAL_NORMAL_X) += data[virialXIndex];
  reduction->item(REDUCTION_VIRIAL_NORMAL_Y) += data[virialYIndex];
  reduction->item(REDUCTION_VIRIAL_NORMAL_Z) += data[virialZIndex];
}


/***************************************************************************
 * RCS INFORMATION:
 *
 *	$RCSfile $
 *	$Author $	$Locker:  $		$State: Exp $
 *	$Revision: 1.1010 $	$Date: 1999/06/17 15:46:01 $
 *
 ***************************************************************************
 * REVISION HISTORY:
 *
 * $Log: ComputeBonds.C,v $
 * Revision 1.1010  1999/06/17 15:46:01  jim
 * Completely rewrote reduction system to eliminate need for sequence numbers.
 *
 * Revision 1.1009  1999/01/06 00:56:20  jim
 * All compute objects except DPMTA now return diagonal of virial tensor.
 *
 * Revision 1.1008  1998/06/18 14:47:59  jim
 * Split virial into NORMAL, NBOND, and SLOW parts to match force classes.
 *
 * Revision 1.1007  1997/10/17 17:16:44  jim
 * Switched from hash tables to checklists, eliminated special exclusion code.
 *
 * Revision 1.1006  1997/09/28 22:36:48  jim
 * Modified tuple-based computations to not duplicate calculations and
 * only require "upstream" proxies.
 *
 * Revision 1.1005  1997/03/20 23:53:29  ari
 * Some changes for comments. Copyright date additions.
 * Hooks for base level update of Compute objects from ComputeMap
 * by ComputeMgr.  Useful for new compute migration functionality.
 *
 * Revision 1.1004  1997/03/19 11:54:03  ari
 * Add Broadcast mechanism.
 * Fixed RCS Log entries on files that did not have Log entries.
 * Added some register variables to Molecule and ComputeNonbondedExcl.C
 *
 ***************************************************************************/
