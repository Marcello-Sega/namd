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
#define MIN_DEBUG_LEVEL 3
#include "Debug.h"

template <class T>
ComputeHomeTuples<T>::ComputeHomeTuples(ComputeID c) : Compute(c) {
  patchMap = PatchMap::Object();
  atomMap = AtomMap::Object();
  reduction = ReductionMgr::Object();

  maxProxyAtoms = 0;
  dummy = NULL;	// initialized to NULL -- won't harm reallocating deletes.
  T::registerReductionData(reduction);
  fake_seq = 0;
}

template <class T>
ComputeHomeTuples<T>::~ComputeHomeTuples()
{
  delete [] dummy;	// allocated during initialize; reallocated in doWork
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
  tuplePatchList.resize(0);

  for ( ai = ai.begin(); ai != ai.end(); ai++ ) {
    tuplePatchList.add(TuplePatchElem((*ai).patch, HOME, cid));
  }

  // Gather all proxy patches (neighbors, that is)
  PatchID neighbors[PatchMap::MaxOneOrTwoAway];
  maxProxyAtoms = 0;
  delete[] dummy;	// deleted, but reallocated soon

  for ( ai = ai.begin(); ai != ai.end(); ai++ ) {
    int numNeighbors = patchMap->oneOrTwoAwayNeighbors((*ai).pid,neighbors);
    for ( int i = 0; i < numNeighbors; ++i )
    {
      if ( patchMap->node(neighbors[i]) != CMyPe() &&
	   ! tuplePatchList.find(TuplePatchElem(neighbors[i])) )
      {
        Patch *patch = patchMap->patch(neighbors[i]);
	tuplePatchList.add(TuplePatchElem(patch, PROXY, cid));
	if (patch->getNumAtoms() > maxProxyAtoms) {
	  maxProxyAtoms = patch->getNumAtoms();
	}
      }
    }
  }

  setNumPatches(tuplePatchList.size());

  // Allocate dummy force vector for non-returning forces
  dummy = new Force[maxProxyAtoms];	// reallocated

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

  tupleList.resize(0);
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
  tupleList.sort(); // This is the expensive operation - tupleList
		    // should be a more efficient container (e.g. hash)
  tupleList.uniq();

  // Resolve all atoms in tupleList to correct PatchList element and index
  ResizeArrayIter<T> al(tupleList);

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
  delete[] dummy;	// deleted, but reallocated very soon
  maxProxyAtoms = 0;

  // find size of largest patch on tuplePatchList, setup dummy force array
  ResizeArrayIter<TuplePatchElem> tpi(tuplePatchList);
  for ( tpi = tpi.begin(); tpi != tpi.end(); tpi++ ) {
    if (tpi->p->getNumAtoms() > maxProxyAtoms) {
      maxProxyAtoms = tpi->p->getNumAtoms();
    }
  }
  dummy = new Force[maxProxyAtoms];	// reallocated
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
  ResizeArrayIter<TuplePatchElem> ap(tuplePatchList);
  for (ap = ap.begin(); ap != ap.end(); ap++) {
    (*ap).x = (*ap).positionBox->open();
    (*ap).a = (*ap).atomBox->open();
    // We only deposit real forces for HomePatch atoms on our Node
    if ( (*ap).patchType == HOME ) 
      (*ap).f = (*ap).forceBox->open();
    else 
      (*ap).f = dummy;
  } 

  BigReal reductionData[T::reductionDataSize];
  for ( int i = 0; i < T::reductionDataSize; ++i ) reductionData[i] = 0;

  // take triplet and pass with tuple info to force eval
  ResizeArrayIter<T> al(tupleList);
  for (al = al.begin(); al != al.end(); al++ ) {
    (*al).computeForce(reductionData);
  }

  T::submitReductionData(reductionData,reduction,fake_seq);
  ++fake_seq;

  // Close boxes - i.e. signal we are done with Positions and
  // AtomProperties and that we are depositing Forces
  for (ap = ap.begin(); ap != ap.end(); ap++) {
    (*ap).positionBox->close(&(*ap).x);
    (*ap).atomBox->close(&(*ap).a);

    if ( (*ap).patchType == HOME ) 
      (*ap).forceBox->close(&(*ap).f);
  }
  DebugM(1, "ComputeHomeTuples::doWork() -- done" << endl);
}
