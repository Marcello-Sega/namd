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

#include "Debug.h"

void NonbondedExclElem::addTuplesForAtom
  (void *voidlist, AtomID atomID, Molecule *molecule)
{
      DebugM(1, "::addTuplesForAtom - atomID " << atomID << endl );
      UniqueSortedArray<NonbondedExclElem> &nonbondedexclList =
                  *( (UniqueSortedArray<NonbondedExclElem>*) voidlist );

      DebugM(1, "::addTuplesForAtom - current list size " << nonbondedexclList.size() << endl );

      /* get list of all nonbondedexcls for the atom */
      LintList *nonbondedexcls = molecule->get_nonbondedexcls_for_atom(atomID);
      DebugM(1, "::addTuplesForAtom - atomID " << atomID << endl );
      DebugM(1, "::addTuplesForAtom - nonbondedexcls->head()" << nonbondedexcls->head() << endl );

      /* cycle through each nonbondedexcl */
      int nonbondedexclNum = nonbondedexcls->head();
      while(nonbondedexclNum != LIST_EMPTY)
      {
        /* store nonbondedexcl in the list */
        DebugM(1, "::addTuplesForAtom - adding nonbondedexcl " << nonbondedexclNum << endl );
        nonbondedexclList.add(NonbondedExclElem(molecule->get_nonbondedexcl(nonbondedexclNum)));
        nonbondedexclNum = nonbondedexcls->next();
      }
}

BigReal NonbondedExclElem::computeForce(void)
{
  DebugM(1, "::computeForce() localIndex = " << localIndex[0] << " "
               << localIndex[1] << endl);

  ComputeNonbondedUtil::calcExcl(
	p[0]->x[localIndex[0]], p[1]->x[localIndex[1]],
	p[0]->f[localIndex[0]], p[1]->f[localIndex[1]],
	p[0]->a[localIndex[0]], p[1]->a[localIndex[1]],
	nonbondedexclType);

  DebugM(3, "::computeForce() -- ending" << endl);
  return(0.);
}

