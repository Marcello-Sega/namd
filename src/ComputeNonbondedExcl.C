/**
***  Copyright (c) 1995, 1996, 1997, 1998, 1999, 2000 by
***  The Board of Trustees of the University of Illinois.
***  All rights reserved.
**/

/*
   Helper for computing non-bonded exclusions 
   Overload of loadTuples() specific for non-bonded exclusions
*/

#include "Namd.h"
#include "Node.h"
#include "Molecule.h"
#include "Parameters.h"
#include "ComputeNonbondedExcl.h"
#include "AtomMap.h"
#include "PatchMap.inl"

//#undef DEBUGM
#include "Debug.h"

void NonbondedExclElem::getMoleculePointers
    (Molecule* mol, int* count, int*** byatom, Exclusion** structarray)
{
  *count = mol->numTotalExclusions;
  *byatom = mol->exclusionsByAtom;
  *structarray = mol->exclusions;
}

void NonbondedExclElem::computeForce(BigReal *reduction)
{
  register TuplePatchElem *p0 = p[0];
  register TuplePatchElem *p1 = p[1];

  register Patch *patch = p0->p;

  if ( patch->flags.doNonbonded )
  {
    register int localIndex0 = localIndex[0];
    register int localIndex1 = localIndex[1];

    nonbonded params;
    params.p_ij = patch->lattice.delta(p0->x[localIndex0], p1->x[localIndex1]);
    params.ff[0] = &(p0->r->f[Results::nbond][localIndex0]);
    params.ff[1] = &(p1->r->f[Results::nbond][localIndex1]);
    params.a[0] = &(p0->a[localIndex0]);
    params.a[1] = &(p1->a[localIndex1]);
    params.m14 = modified;
    params.reduction = reduction;

    if ( patch->flags.doFullElectrostatics )
    {
      params.fullf[0] = &(p0->r->f[Results::slow][localIndex0]);
      params.fullf[1] = &(p1->r->f[Results::slow][localIndex1]);
      if ( patch->flags.doMolly ) {
        ComputeNonbondedUtil::calcExcl(&params);
        params.p_ij =
	  patch->lattice.delta(p0->x_avg[localIndex0], p1->x_avg[localIndex1]);
        ComputeNonbondedUtil::calcSlowExcl(&params);
      } else {
        ComputeNonbondedUtil::calcFullExcl(&params);
      }
    }
    else
      ComputeNonbondedUtil::calcExcl(&params);
  }
}

#if(0)

#include "charm++.h"
#include "Inform.h"
#include "Node.h"
#include "Molecule.h"
#include "Parameters.h"
#include "Node.h"
#include "UniqueSet.h"
#include "UniqueSetIter.h"

void
ComputeNonbondedExcls::loadTuples() {

  // cycle through each patch and gather all tuples
  TuplePatchListIter ai(tuplePatchList);
  register int i;
  int numExclusions = node->molecule->numTotalExclusions;

  char *exclFlag = new char[numExclusions];
  for (register char *c = exclFlag; c < (exclFlag+numExclusions); *c++ = 0);

  tupleList.clear();
  for ( ai = ai.begin(); ai != ai.end(); ai++ )
  {
    Patch *patch = (*ai).p;
    AtomIDList atomID = patch->getAtomIDList();

    // cycle through each atom in the patch and load up tuples
    int numAtoms = patch->getNumAtoms();
    for (i=0; i < numAtoms; i++)
    {
       /* get list of all bonds for the atom */
       register int *excls = node->molecule->get_exclusions_for_atom(atomID[i]);

       /* cycle through each exclusion */
       while(*excls != -1) {
         exclFlag[*excls++] = 1;
       }
    }
  }

  if ( node->simParameters->fixedAtomsOn ) {
    Molecule *molecule = node->molecule;
    for (i=0; i<numExclusions; i++) {
      if (exclFlag[i]) {
	Exclusion *excl = molecule->get_exclusion(i);
	if ( ! ( molecule->is_atom_fixed(excl->atom1) &&
		 molecule->is_atom_fixed(excl->atom2) ) )
	  tupleList.load(NonbondedExclElem(excl));
      }
    }
  } else {
    for (i=0; i<numExclusions; i++) {
      if (exclFlag[i]) {
	tupleList.load(NonbondedExclElem(node->molecule->get_exclusion(i)));
      }
    }
  }
  delete[] exclFlag;
 
  // Resolve all atoms in tupleList to correct PatchList element and index
  // and eliminate tuples we aren't responsible for
  UniqueSetIter<NonbondedExclElem> al(tupleList);
  UniqueSet<NonbondedExclElem> tupleList2;

  LocalID aid[NonbondedExclElem::size];
  for (al = al.begin(); al != al.end(); al++ ) {
    aid[0] = atomMap->localID(al->atomID[0]);
    aid[1] = atomMap->localID(al->atomID[1]);
    int homepatch = patchMap->downstream(aid[0].pid,aid[1].pid);;
    NonbondedExclElem &t = *al;
    if ( homepatch != notUsed && patchMap->node(homepatch) == CkMyPe() ) {
      t.p[0] = tuplePatchList.find(TuplePatchElem(aid[0].pid));
      t.p[1] = tuplePatchList.find(TuplePatchElem(aid[1].pid));
      t.localIndex[0] = aid[0].index;
      t.localIndex[1] = aid[1].index;
      tupleList2.load(t);
    }
  }
  tupleList = tupleList2;

}

#endif

