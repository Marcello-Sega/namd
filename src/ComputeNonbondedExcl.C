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

#include "ComputeNonbondedExcl.h"
#include "Molecule.h"
#include "Parameters.h"
#include "Node.h"

#undef DEBUGM
#include "Debug.h"

void NonbondedExclElem::loadTuplesForAtom
  (void *voidlist, AtomID atomID, Molecule *molecule)
{
      DebugM(1, "::loadTuplesForAtom - atomID " << atomID << endl );
      UniqueSortedArray<NonbondedExclElem> &exclList =
                  *( (UniqueSortedArray<NonbondedExclElem>*) voidlist );

      DebugM(1, "::loadTuplesForAtom - current list size " << exclList.size() << endl );

      /* get list of all bonds for the atom */
      LintList *excls = molecule->get_exclusions_for_atom(atomID);
      DebugM(1, "::loadTuplesForAtom - atomID " << atomID << endl );
      DebugM(1, "::loadTuplesForAtom - excls->head()" << excls->head() << endl );

      /* cycle through each exclusion */
      int exclNum = excls->head();
      while(exclNum != LIST_EMPTY)
      {
        /* store exclusion in the list */
        DebugM(1, "::loadTuplesForAtom - adding excl " << exclNum << endl );
        exclList.load(NonbondedExclElem(molecule->get_exclusion(exclNum)));
        exclNum = excls->next();
      }
}

BigReal NonbondedExclElem::reductionDummy[reductionDataSize];

void NonbondedExclElem::computeForce(BigReal *reduction)
{
  register TuplePatchElem *p0 = p[0];
  register TuplePatchElem *p1 = p[1];

  if ( p0->patchType != HOME ) reduction = reductionDummy;

  register Patch *patch = p0->p;

  register int localIndex0 = localIndex[0];
  register int localIndex1 = localIndex[1];

  Vector x01(patch->lattice.delta(p0->x[localIndex0], p1->x[localIndex1]));

  if ( patch->flags.doFullElectrostatics )
    ComputeNonbondedUtil::calcFullExcl(
	x01,
	p0->f[localIndex0], p1->f[localIndex1],
	p0->f[localIndex0], p1->f[localIndex1],
	p0->a[localIndex0], p1->a[localIndex1],
	modified, reduction);
  else
    ComputeNonbondedUtil::calcExcl(
	x01,
	p0->f[localIndex0], p1->f[localIndex1],
	p0->a[localIndex0], p1->a[localIndex1],
	modified, reduction);
}


