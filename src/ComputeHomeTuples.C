/***************************************************************************/
/*           (C) Copyright 1996,1997 The Board of Trustees of the          */
/*                          University of Illinois                         */
/*                           All Rights Reserved                           */
/***************************************************************************/
/***************************************************************************
 * DESCRIPTION: General template class for computing Bonded Forces 
 *              like Angles, Dihedrals etc..
 *
 ***************************************************************************/

#include "Namd.h"
#include "Node.h"
#include "PatchMap.h"
#include "AtomMap.h"
#include "ComputeHomeTuples.h"
#include "PatchMgr.h"
#include "HomePatchList.h"
#include "Molecule.h"
#include "ReductionMgr.h"
#include "Inform.h"
#include "Templates/UniqueSet.h"
#include "Templates/UniqueSetIter.h"

#define DEBUGM
#undef MIN_DEBUG_LEVEL
#define MIN_DEBUG_LEVEL 3
#include "Debug.h"

template <class T>
ComputeHomeTuples<T>::ComputeHomeTuples(ComputeID c) : Compute(c) {
  patchMap = PatchMap::Object();
  atomMap = AtomMap::Object();
  reduction = ReductionMgr::Object();

  maxProxyAtoms = 0;
  dummyForce = NULL;	// initialized to NULL -- won't harm reallocating deletes.
  T::registerReductionData(reduction);
}

template <class T>
ComputeHomeTuples<T>::~ComputeHomeTuples()
{
  delete [] dummyForce;	// allocated during initialize; reallocated in doWork
  T::unregisterReductionData(reduction);
}


//===========================================================================
// initialize() - Method is invoked only the first time
// atom maps, patchmaps etc are ready and we are about to start computations
//===========================================================================
template <class T>
void ComputeHomeTuples<T>::initialize() {

  // Gather all HomePatches
  HomePatchList *a = patchMap->homePatchList();
  ResizeArrayIter<HomePatchElem> ai(*a);

  // Start with empty list
  tuplePatchList.clear();

  for ( ai = ai.begin(); ai != ai.end(); ai++ ) {
    tuplePatchList.add(TuplePatchElem((*ai).patch, HOME, cid));
  }

  // Gather all proxy patches (neighbors, that is)
  PatchID neighbors[PatchMap::MaxOneOrTwoAway];

  for ( ai = ai.begin(); ai != ai.end(); ai++ ) {
    int numNeighbors = patchMap->oneOrTwoAwayNeighbors((*ai).pid,neighbors);
    for ( int i = 0; i < numNeighbors; ++i )
    {
      if ( patchMap->node(neighbors[i]) != CMyPe() &&
	   ! tuplePatchList.find(TuplePatchElem(neighbors[i])) )
      {
        Patch *patch = patchMap->patch(neighbors[i]);
	tuplePatchList.add(TuplePatchElem(patch, PROXY, cid));
      }
    }
  }

  setNumPatches(tuplePatchList.size());

  sizeDummy();
  loadTuples();
}

//===========================================================================
// atomUpdate() - Method is invoked after anytime that atoms have been
// changed in patches used by this Compute object.
//===========================================================================
template <class T>
void ComputeHomeTuples<T>::atomUpdate() {
  sizeDummy();
  loadTuples();
}

template <class T>
void ComputeHomeTuples<T>::loadTuples() {
  // cycle through each home patch and gather all tuples
  HomePatchList *a = patchMap->homePatchList();
  ResizeArrayIter<HomePatchElem> ai(*a);

  tupleList.clear();
  for ( ai = ai.begin(); ai != ai.end(); ai++ )
  {
    Patch *patch = (*ai).patch;
    AtomIDList atomID = patch->getAtomIDList();

    // cycle through each atom in the patch and load up tuples
    for (int i=0; i < patch->getNumAtoms(); i++)
    {
      T::loadTuplesForAtom((void*)&tupleList,atomID[i],node->molecule);
    }
  }
  tupleList.rehash();

  // Resolve all atoms in tupleList to correct PatchList element and index
  UniqueSetIter<T> al(tupleList);

  for (al = al.begin(); al != al.end(); al++ ) {
    for (int i=0; i < T::size; i++) {
	LocalID aid = atomMap->localID(al->atomID[i]);
	al->p[i] = tuplePatchList.find(TuplePatchElem(aid.pid));
	if ( ! (al->p)[i] )
	{
	  iout << iERROR << "ComputeHomeTuples couldn't find patch " 
	    << aid.pid << " for atom " << al->atomID[i] 
	    << ", aborting.\n" << endi;
	  Namd::die();
	}
	al->localIndex[i] = aid.index;
    }
  }
}


//-----------------------------------------------------------------------
// Figure out maximum # of atoms in the Proxies we will be dealing with
// and allocate a dummy force vector of that size to take in
// the force calculation (which is thrown away since only HomePatch
// forces are of interest)
//-----------------------------------------------------------------------
template <class T>
void ComputeHomeTuples<T>::sizeDummy() {
  delete[] dummyForce;	// deleted, but reallocated very soon
  maxProxyAtoms = 0;

  // find size of largest patch on tuplePatchList, setup dummy force array
  UniqueSetIter<TuplePatchElem> tpi(tuplePatchList);
  for ( tpi = tpi.begin(); tpi != tpi.end(); tpi++ ) {
    if (tpi->p->getNumAtoms() > maxProxyAtoms) {
      maxProxyAtoms = tpi->p->getNumAtoms();
    }
  }
  dummyForce = new Force[maxProxyAtoms];	// reallocated
  for ( int i = 0; i < Results::maxNumForces; ++i )
  {
    dummyResults.f[i] = dummyForce;
  }
}


//-------------------------------------------------------------------
// Routine which is called by enqueued work msg.  It wraps
// actualy Force computation with the apparatus needed
// to get access to atom positions, return forces etc.
//-------------------------------------------------------------------
template <class T>
void ComputeHomeTuples<T>::doWork() {
  DebugM(1, "ComputeHomeTuples::doWork() -- started " << endl );

  // Open Boxes - register tFat we are using Positions
  // and will be depositing Forces.
  UniqueSetIter<TuplePatchElem> ap(tuplePatchList);
  for (ap = ap.begin(); ap != ap.end(); ap++) {
    ap->x = ap->positionBox->open();
    ap->a = ap->atomBox->open();
    // We only deposit real forces for HomePatch atoms on our Node
    if ( ap->patchType == HOME ) 
    {
      ap->r = ap->forceBox->open();
    }
    else 
    {
      ap->r = &dummyResults;
    }
    ap->f = ap->r->f[Results::normal];
  } 

  BigReal reductionData[T::reductionDataSize];
  for ( int i = 0; i < T::reductionDataSize; ++i ) reductionData[i] = 0;

  // take triplet and pass with tuple info to force eval
  UniqueSetIter<T> al(tupleList);
  for (al = al.begin(); al != al.end(); al++ ) {
    al->computeForce(reductionData);
  }

  T::submitReductionData(reductionData,reduction,ap.begin()->p->flags.seq);

  // Close boxes - i.e. signal we are done with Positions and
  // AtomProperties and that we are depositing Forces
  for (ap = ap.begin(); ap != ap.end(); ap++) {
    ap->positionBox->close(&(ap->x));
    ap->atomBox->close(&(ap->a));

    if ( ap->patchType == HOME ) 
      ap->forceBox->close(&(ap->r));
  }
  DebugM(1, "ComputeHomeTuples::doWork() -- done" << endl);
}


/***************************************************************************
 * RCS INFORMATION:
 *
 *      $RCSfile: ComputeHomeTuples.C,v $
 *      $Author: milind $  $Locker:  $             $State: Exp $
 *      $Revision: 1.1010 $     $Date: 1997/04/04 23:34:16 $
 *
 ***************************************************************************
 * REVISION HISTORY:
 *
 * $Log: ComputeHomeTuples.C,v $
 * Revision 1.1010  1997/04/04 23:34:16  milind
 * Got NAMD2 to run on Origin2000.
 * Included definitions of class static variables in C files.
 * Fixed alignment bugs by using memcpy instead of assignment in
 * pack and unpack.
 *
 * Revision 1.1009  1997/03/18 21:35:25  jim
 * Eliminated fake_seq.  Reductions now use Patch::flags.seq.
 *
 * Revision 1.1008  1997/03/13 22:39:35  jim
 * Fixed some bugs in multiple-force return / full electrostatics.
 *
 * Revision 1.1007  1997/03/12 22:06:35  jim
 * First step towards multiple force returns and multiple time stepping.
 *
 * Revision 1.1006  1997/03/11 23:46:26  ari
 * Improved ComputeNonbondedExcl loadTuples() by overloading the default
 * template method from ComputeHomeTuples and used the checklist suggested
 * by Jim.  Good performance gain.
 *
 * Revision 1.1005  1997/03/10 17:40:05  ari
 * UniqueSet changes - some more commenting and cleanup
 *
 * Revision 1.1004  1997/03/04 22:38:16  ari
 * Reworked ResizeArray - much more rational.  Overall about
 * same performance as before (sometimes a little better).
 * Needed tricks to make ResizeArray(Raw) work fast.
 * Clean up of code.
 *
 *
 ***************************************************************************/
