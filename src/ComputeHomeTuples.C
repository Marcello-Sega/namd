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
  DebugM(1, "ComputeHomeTuples::ComputeHomeTuples(%d) -- starting " << (int)c << endl );
  patchMap = PatchMap::Object();
  atomMap = AtomMap::Object();
  reduction = ReductionMgr::Object();

  maxProxyAtoms = 0;
  dummy = NULL;
  T::registerReductionData(reduction);
  fake_seq = 0;

  DebugM(1, "ComputeHomeTuples::ComputeHomeTuples(%d) -- done " << (int)c << endl);
}

template <class T>
ComputeHomeTuples<T>::~ComputeHomeTuples()
{
  delete [] dummy;
  T::unregisterReductionData(reduction);
}


//===========================================================================
// initialize() - Method is invoked only the first time
// atom maps, patchmaps etc are ready and we are about to start computations
//===========================================================================
template <class T>
void ComputeHomeTuples<T>::initialize() {

  // Gather all patches
  DebugM(4, "ComputeHomeTuples::initialize() - Starting Up" << endl );
  HomePatchList *a = patchMap->homePatchList();
  ResizeArrayIter<HomePatchElem> ai(*a);

  tuplePatchList.resize(0);
  DebugM(1, "ComputeHomeTuples::initialize() - Size of the tuplePatchList " << tuplePatchList.size() << endl );

  for ( ai = ai.begin(); ai != ai.end(); ai++ ) {
    tuplePatchList.add(TuplePatchElem((*ai).p, HOME, cid));
    DebugM( 1, "ComputeHomeTuples::initialize() - adding Patch " << (*ai).p->getPatchID() << " to list" << endl );
  }

  // Gather all proxy patches (neighbors, that is)
  PatchID neighbors[PatchMap::MaxOneOrTwoAway];
  maxProxyAtoms = 0;
  delete[] dummy;

  for ( ai = ai.begin(); ai != ai.end(); ai++ ) {
    int numNeighbors = patchMap->oneOrTwoAwayNeighbors((*ai).pid,neighbors);
    for ( int i = 0; i < numNeighbors; ++i )
    {
      if ( patchMap->node(neighbors[i]) != CMyPe() &&
	   ! tuplePatchList.find(TuplePatchElem(neighbors[i])) )
      {
        Patch *patch = patchMap->patch(neighbors[i]);
	DebugM( 1, "ComputeHomeTuples::initialize() - adding (Proxy)Patch " <<
		patch->getPatchID() << " to list" << endl );
	tuplePatchList.add(TuplePatchElem(patch, PROXY, cid));
	DebugM( 1, "ComputeHomeTuples::initialize() - tuplePatchList now has " <<
		tuplePatchList.size() << " elements" << endl );
	if (patch->getNumAtoms() > maxProxyAtoms) {
	  maxProxyAtoms = patch->getNumAtoms();
	}
      }
    }
  }

  setNumPatches(tuplePatchList.size());

  dummy = new Force[maxProxyAtoms];

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
    Patch *p = (*ai).p;
    AtomIDList atomID = p->getAtomIDList();

    // cycle through each atom in the patch and load up tuples
    for (int i=0; i < p->getNumAtoms(); i++)
    {
      T::loadTuplesForAtom((void*)&tupleList,atomID[i],node->molecule);
    }
  }
  tupleList.sort();
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


template <class T>
void ComputeHomeTuples<T>::sizeDummy() {
  delete[] dummy;
  maxProxyAtoms = 0;

  // find size of largest patch on tuplePatchList, setup dummy force array
  ResizeArrayIter<TuplePatchElem> tpi;
  for ( tpi = tpi.begin(); tpi != tpi.end(); tpi++ ) {
    if (tpi->p->getNumAtoms() > maxProxyAtoms) {
      maxProxyAtoms = tpi->p->getNumAtoms();
    }
  }
  dummy = new Force[maxProxyAtoms];
}


template <class T>
void ComputeHomeTuples<T>::doWork() {
  DebugM(1, "ComputeHomeTuples::doWork() -- started " << endl );
  // Open Boxes
  ResizeArrayIter<TuplePatchElem> ap(tuplePatchList);
  for (ap = ap.begin(); ap != ap.end(); ap++) {
    (*ap).x = (*ap).positionBox->open();
    if ( (*ap).patchType == HOME ) (*ap).f = (*ap).forceBox->open();
    else (*ap).f = dummy;
    (*ap).a = (*ap).atomBox->open();
  } 

  BigReal reductionData[T::reductionDataSize];
  for ( int i = 0; i < T::reductionDataSize; ++i ) reductionData[i] = 0;

  // take triplet and pass with tuple info to force eval
  DebugM(3, "ComputeHomeTuples::doWork() - size of tuple list = " << tupleList.size() << endl );
  ResizeArrayIter<T> al(tupleList);
  for (al = al.begin(); al != al.end(); al++ ) {
    (*al).computeForce(reductionData);
  }

  T::submitReductionData(reductionData,reduction,fake_seq);
  ++fake_seq;

  for (ap = ap.begin(); ap != ap.end(); ap++) {
    (*ap).positionBox->close(&(*ap).x);
    if ( (*ap).patchType == HOME ) (*ap).forceBox->close(&(*ap).f);
    (*ap).atomBox->close(&(*ap).a);
  }
  DebugM(1, "ComputeHomeTuples::doWork() -- done" << endl);
}

