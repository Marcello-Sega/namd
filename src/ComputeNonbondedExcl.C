/**
***  Copyright (c) 1995, 1996, 1997, 1998, 1999, 2000 by
***  The Board of Trustees of the University of Illinois.
***  All rights reserved.
**/

/*
   Helper for computing non-bonded exclusions 
   Overload of loadTuples() specific for non-bonded exclusions
*/

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

