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

void NonbondedExclElem::addTuplesForAtom
  (void *voidlist, AtomID atomID, Molecule *molecule)
{
      DebugM(1, "::addTuplesForAtom - atomID " << atomID << endl );
      UniqueSortedArray<NonbondedExclElem> &exclList =
                  *( (UniqueSortedArray<NonbondedExclElem>*) voidlist );

      DebugM(1, "::addTuplesForAtom - current list size " << exclList.size() << endl );

      /* get list of all bonds for the atom */
      LintList *excls = molecule->get_exclusions_for_atom(atomID);
      DebugM(1, "::addTuplesForAtom - atomID " << atomID << endl );
      DebugM(1, "::addTuplesForAtom - excls->head()" << excls->head() << endl );

      /* cycle through each exclusion */
      int exclNum = excls->head();
      while(exclNum != LIST_EMPTY)
      {
        /* store exclusion in the list */
        DebugM(1, "::addTuplesForAtom - adding excl " << exclNum << endl );
        exclList.add(NonbondedExclElem(molecule->get_exclusion(exclNum)));
        exclNum = excls->next();
      }
}


void NonbondedExclElem::computeForce(BigReal *reduction)
{
  DebugM(1, "::computeForce() localIndex = " << localIndex[0] << " "
               << localIndex[1] << endl);

  BigReal dummy[reductionDataSize];

  ComputeNonbondedUtil::calcExcl(
	p[0]->x[localIndex[0]], p[1]->x[localIndex[1]],
	p[0]->f[localIndex[0]], p[1]->f[localIndex[1]],
	p[0]->a[localIndex[0]], p[1]->a[localIndex[1]],
	modified,
	( p[0]->patchType == HOME ) ? reduction : dummy );

  DebugM(3, "::computeForce() -- ending" << endl);
}


