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
#include "ComputeNonbondedUtil.h"
#include "Molecule.h"
#include "Parameters.h"
#include "Node.h"
#define DEBUGM
#include "Debug.h"

void NonbondedExclElem::addTuplesForAtom
  (void *voidlist, AtomID atomID, Molecule *molecule)
{
      DebugM(1, "::addTuplesForAtom - atomID " << atomID << endl );
      UniqueSortedArray<NonbondedExclElem> &nonbondedexclList =
                  *( (UniqueSortedArray<NonbondedExclElem>*) voidlist );
      NonbondedExclElem E;

      DebugM(1, "::addTuplesForAtom - current list size " << nonbondedexclList.size() << endl );

      /* cycle through each regular exclusion */
      IntList *nonbondedexcls =molecule->get_nonbondedexcls_for_allatom(atomID);
      DebugM(1, "::addTuplesForAtom - atomID " << atomID << endl );
      int maxNum = nonbondedexcls->num();
      int nonbondedexclNum;
      for(nonbondedexclNum = 0; nonbondedexclNum < maxNum; nonbondedexclNum++)
      {
        /* store nonbondedexcl in the list */
        DebugM(1,"::addTuplesForAtom - adding nonbondedexcl "
		 << nonbondedexclNum << " (regular type)" << "\n" );
DebugM(2,"Adding " << NonbondedExclElem(molecule->get_nonbondedexcl(nonbondedexclNum)).atomID[0] << " " << NonbondedExclElem(molecule->get_nonbondedexcl(nonbondedexclNum)).atomID[1]<< "\n");
	E = NonbondedExclElem(molecule->get_nonbondedexcl(nonbondedexclNum));
DebugM(2,"Adding " << E.atomID[0] << " " << E.atomID[1] << "\n");
	E.nonbondedexclType = 0;
        nonbondedexclList.add(E);
      }

      DebugM(1,"::addTuplesForAtom - almost done size is "
		<< nonbondedexclList.size() << "\n");

      /* cycle through each (special) 14 exclusion */
      nonbondedexcls = molecule->get_nonbondedexcls_for_14atom(atomID);
      DebugM(1, "::addTuplesForAtom - atomID " << atomID << "\n" );
      maxNum = nonbondedexcls->num();
      for(nonbondedexclNum = 0; nonbondedexclNum < maxNum; nonbondedexclNum++)
      {
        /* store nonbondedexcl in the list */
        DebugM(1,"::addTuplesForAtom - adding nonbondedexcl "
		 << nonbondedexclNum << " (onefour type)" << "\n" );
DebugM(2,"Adding " << NonbondedExclElem(molecule->get_nonbondedexcl(nonbondedexclNum)).atomID[0] << " " << NonbondedExclElem(molecule->get_nonbondedexcl(nonbondedexclNum)).atomID[1]<< "\n");
	E = NonbondedExclElem(molecule->get_nonbondedexcl(nonbondedexclNum));
	E.nonbondedexclType = 1;
        nonbondedexclList.add(E);
      }
      DebugM(1,"::addTuplesForAtom - done size is "
		<< nonbondedexclList.size() << "\n");
}

BigReal NonbondedExclElem::computeForce(void)
{
  DebugM(1, "::computeForce() localIndex = " << localIndex[0] << " "
               << localIndex[1] << " type: " << nonbondedexclType  << "\n");

  ComputeNonbondedUtil::calcExcl(
	p[0]->x[localIndex[0]], p[1]->x[localIndex[1]],
	p[0]->f[localIndex[0]], p[1]->f[localIndex[1]],
	p[0]->a[localIndex[0]], p[1]->a[localIndex[1]],
	nonbondedexclType);

  DebugM(3, "::computeForce() -- ending" << endl);
  return(0.);
}

