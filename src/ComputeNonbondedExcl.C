/***************************************************************************/
/*         (C) Copyright 1996,1997 The Board of Trustees of the            */
/*                          University of Illinois                         */
/*                           All Rights Reserved                           */
/***************************************************************************/
/***************************************************************************
 * DESCRIPTION: Helper for computing non-bonded exclusions 
 *		Overload of loadTuples() specific for non-bonded exclusions
 *
 ***************************************************************************/

#include "Namd.h"
#include "Node.h"
#include "Molecule.h"
#include "Parameters.h"
#include "ComputeNonbondedExcl.h"
#include "AtomMap.h"

//#undef DEBUGM
#include "Debug.h"

BigReal NonbondedExclElem::reductionDummy[reductionDataSize];

void NonbondedExclElem::computeForce(BigReal *reduction)
{
  register TuplePatchElem *p0 = p[0];
  register TuplePatchElem *p1 = p[1];

  if ( p0->patchType != HOME ) reduction = reductionDummy;

  register Patch *patch = p0->p;

  if ( patch->flags.doNonbonded )
  {
    register int localIndex0 = localIndex[0];
    register int localIndex1 = localIndex[1];

    Vector x01(patch->lattice.delta(p0->x[localIndex0], p1->x[localIndex1]));

    if ( patch->flags.doFullElectrostatics )
      ComputeNonbondedUtil::calcFullExcl(
	x01,
	p0->r->f[Results::nbond][localIndex0],
	p1->r->f[Results::nbond][localIndex1],
	p0->r->f[Results::slow][localIndex0],
	p1->r->f[Results::slow][localIndex1],
	p0->a[localIndex0], p1->a[localIndex1],
	modified, reduction);
    else
      ComputeNonbondedUtil::calcExcl(
	x01,
	p0->r->f[Results::nbond][localIndex0],
	p1->r->f[Results::nbond][localIndex1],
	p0->a[localIndex0], p1->a[localIndex1],
	modified, reduction);
  }
}


#include "ckdefs.h"
#include "chare.h"
#include "c++interface.h"
#include "Inform.h"
#include "Node.h"
#include "Molecule.h"
#include "Parameters.h"
#include "Node.h"
#include "Templates/UniqueSet.h"
#include "Templates/UniqueSetIter.h"

void
ComputeNonbondedExcls::loadTuples() {

  // cycle through each home patch and gather all tuples
  HomePatchList *a = patchMap->homePatchList();
  ResizeArrayIter<HomePatchElem> ai(*a);
  register int i;
  int numExclusions = node->molecule->numTotalExclusions;

  char *exclFlag = new char[numExclusions];
  for (register char *c = exclFlag; c < (exclFlag+numExclusions); *c++ = 0);

  tupleList.clear();
  for ( ai = ai.begin(); ai != ai.end(); ai++ )
  {
    Patch *patch = (*ai).patch;
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
  for (i=0; i<numExclusions; i++) {
    if (exclFlag[i]) {
      tupleList.load(NonbondedExclElem(node->molecule->get_exclusion(i)));
    }
  }
  delete[] exclFlag;
 
  // Resolve all atoms in tupleList to correct PatchList element and index
  UniqueSetIter<NonbondedExclElem> al(tupleList);
 
  for (al = al.begin(); al != al.end(); al++ ) {
    for (i=0; i < NonbondedExclElem::size; i++) {
      LocalID aid = atomMap->localID(al->atomID[i]);
      al->p[i] = tuplePatchList.find(TuplePatchElem(aid.pid));
      /*
      if ( ! (al->p)[i] ) {
 	iout << iERROR << "ComputeHomeTuples couldn't find patch " 
 	    << aid.pid << " for atom " << al->atomID[i] 
 	    << ", aborting.\n" << endi;
 	Namd::die();
      }
      */
      al->localIndex[i] = aid.index;
    }
  }
}

/***************************************************************************
 * RCS INFORMATION:
 *
 *	$RCSfile: ComputeNonbondedExcl.C,v $
 *	$Author: ari $	$Locker:  $		$State: Exp $
 *	$Revision: 1.1012 $	$Date: 1997/04/06 22:44:59 $
 *
 ***************************************************************************
 * REVISION HISTORY:
 *
 * $Log: ComputeNonbondedExcl.C,v $
 * Revision 1.1012  1997/04/06 22:44:59  ari
 * Add priorities to messages.  Mods to help proxies without computes.
 * Added quick enhancement to end of list insertion of ResizeArray(s)
 *
 * Revision 1.1011  1997/03/25 23:00:56  jim
 * Added nonbondedFrequency parameter and multiple time-stepping
 *
 * Revision 1.1010  1997/03/19 11:54:11  ari
 * Add Broadcast mechanism.
 * Fixed RCS Log entries on files that did not have Log entries.
 * Added some register variables to Molecule and ComputeNonbondedExcl.C
 *
 * Revision 1.1009  1997/03/19 05:49:59  jim
 * Added ComputeSphericalBC, cleaned up make dependencies.
 *
 * Revision 1.1008  1997/03/13 06:36:57  jim
 * Multiple time-stepping implemented, still needs proper splitting functions.
 *
 * Revision 1.1007  1997/03/11 23:46:28  ari
 * Improved ComputeNonbondedExcl loadTuples() by overloading the default
 * template method from ComputeHomeTuples and used the checklist suggested
 * by Jim.  Good performance gain.
 *
 * Revision 1.1002  1997/03/10 17:40:09  ari
 * UniqueSet changes - some more commenting and cleanup
 *
 * Revision 1.1001  1997/03/09 22:28:24  jim
 * (Hopefully) sped up exclusion calculation (removed gross inefficiencies).
 *
 * Revision 1.1000  1997/02/06 15:58:08  ari
 * Resetting CVS to merge branches back into the main trunk.
 * We will stick to main trunk development as suggested by CVS manual.
 * We will set up tags to track fixed points of development/release
 * as suggested by CVS manual - all praise the CVS manual.
 *
 * Revision 1.779  1997/02/06 15:53:04  ari
 * Updating Revision Line, getting rid of branches
 *
 * Revision 1.778.2.1  1997/02/05 22:18:08  ari
 * Added migration code - Currently the framework is
 * there with compiling code.  This version does
 * crash shortly after migration is complete.
 * Migration appears to complete, but Patches do
 * not appear to be left in a correct state.
 *
 * Revision 1.778  1997/01/28 00:30:20  ari
 * internal release uplevel to 1.778
 *
 * Revision 1.777.2.1  1997/01/27 22:45:05  ari
 * Basic Atom Migration Code added.
 * Added correct magic first line to .h files for xemacs to go to C++ mode.
 * Compiles and runs without migration turned on.
 *
 * Revision 1.777  1997/01/17 19:35:54  ari
 * Internal CVS leveling release.  Start development code work
 * at 1.777.1.1.
 *
 * Revision 1.4  1997/01/16 20:00:03  jim
 * Added reduction calls to ComputeNonbondedSelf and ...Pair.
 * Also moved some code from ...Excl to ...Util.
 *
 * Revision 1.3  1997/01/16 00:56:01  jim
 * Added reduction of energies from ComputeHomeTuples objects, except
 * for ComputeNonbondedExcl which only reports 0 energy.
 * Some problems with ReductionMgr are apparent, but it still runs.
 *
 * Revision 1.2  1996/12/06 06:56:11  jim
 * cleaned up and renamed a bit, now it works
 *
 * Revision 1.1  1996/12/03 17:17:03  nealk
 * Initial revision
 *
 * Revision 1.2  1996/12/03 15:15:40  nealk
 * Removed tons-o-debugging.
 *
 * Revision 1.1  1996/12/03 14:53:42  nealk
 * Initial revision
 *
 *
 ***************************************************************************/

