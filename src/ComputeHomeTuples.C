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

#include "Node.h"
#include "PatchMap.h"
#include "AtomMap.h"
#include "ComputeHomeTuples.h"
#include "PatchMgr.h"
#include "Molecule.h"
#define DEBUGM
#define MIN_DEBUG_LEVEL 1
#include "Debug.h"

template <class T>
ComputeHomeTuples<T>::ComputeHomeTuples(ComputeID c) : Compute(c) {
  DebugM(1, "ComputeHomeTuples::ComputeHomeTuples(%d) -- starting " << (int)c << endl );
  patchMap = PatchMap::Object();
  atomMap = AtomMap::Object();

  maxProxyAtoms = 0;
  dummy = NULL;
  DebugM(1, "ComputeHomeTuples::ComputeHomeTuples(%d) -- done " << (int)c << endl);
}

template <class T>
void ComputeHomeTuples<T>::mapReady() {

  // Gather all patches
  DebugM(1, "ComputeHomeTuples::mapReady() - Starting Up" << endl );
  HomePatchList *a = patchMap->homePatchList();
  ResizeArrayIter<HomePatchElem> ai(*a);

  setNumPatches(a->size());
  tuplePatchList.resize(0);
  DebugM(1, "ComputeHomeTuples::mapReady() - Size of the tuplePatchList " << tuplePatchList.size() << endl );

  for ( ai = ai.begin(); ai != ai.end(); ai++ ) {
    tuplePatchList.add(TuplePatchElem((*ai).p, HOME, cid));
    DebugM( 1, "ComputeHomeTuples::mapReady() - adding Patch " << (*ai).p->getPatchID() << " to list" << endl );
  }

  // Gather all proxy patches (neighbors, that is)
  PatchID oneaway[PatchMap::MaxOneAway];
  maxProxyAtoms = 0;
  delete[] dummy;

  for ( ai = ai.begin(); ai != ai.end(); ai++ ) {
    int numOneAway = patchMap->oneAwayNeighbors((*ai).pid, oneaway);
    for ( int i = 0; i < numOneAway; ++i )
    {
      if ( patchMap->node(oneaway[i]) != CMyPe() &&
	   ! tuplePatchList.find(TuplePatchElem(oneaway[i])) )
      {
        Patch *patch = patchMap->patch(oneaway[i]);
	DebugM( 1, "ComputeHomeTuples::mapReady() - adding (Proxy)Patch " <<
		patch->getPatchID() << " to list" << endl );
	tuplePatchList.add(TuplePatchElem(patch, PROXY, cid));
	DebugM( 1, "ComputeHomeTuples::mapReady() - tuplePatchList now has " <<
		tuplePatchList.size() << " elements" << endl );
	if (patch->getNumAtoms() > maxProxyAtoms) {
	  maxProxyAtoms = patch->getNumAtoms();
	}
      }
    }
  }

  dummy = new Force[maxProxyAtoms];

  ResizeArrayIter<TuplePatchElem> tpi(tuplePatchList);

  /* cycle through each patch */
  DebugM(1, "ComputeHomeTuples::mapReady() - iterating over patches to get atoms" << endl);
  for ( tpi = tpi.begin(); tpi != tpi.end(); tpi++ )
  {
    Patch *p = (*tpi).p;
    DebugM(1, "ComputeHomeTuples::mapReady() - looking at patch " <<
	p->getPatchID() << " with " << p->getNumAtoms() << " atoms" << endl );
    AtomIDList atomID = p->getAtomIDList();
    DebugM(1, "ComputeHomeTuples::mapReady() - confirm patch " <<
	p->getPatchID() << " with " << atomID.size() << " atoms" << endl );

    /* cycle through each atom in the patch */

    for (int i=0; i < p->getNumAtoms(); i++)
    {
      T::addTuplesForAtom((void*)&tupleList,atomID[i],node->molecule);
    }
    DebugM(1, "ComputeHomeTuples::mapReady() - tupleList size = " << tupleList.size() << endl );
  }

  // Resolve all atoms in tupleList to correct PatchList element and index
  ResizeArrayIter<T> al(tupleList);

  for (al = al.begin(); al != al.end(); al++ ) {
    for (int i=0; i < T::size; i++) {
	LocalID aid = atomMap->localID((*al).atomID[i]);
	(*al).p[i] = tuplePatchList.find(TuplePatchElem(aid.pid));
	(*al).localIndex[i] = aid.index;
    }
  }
}


template <class T>
void ComputeHomeTuples<T>::doWork() {
  DebugM(1, "ComputeHomeTuples::doWork() -- started " << endl );
  // Open Boxes
  ResizeArrayIter<TuplePatchElem> ap(tuplePatchList);
  for (ap = ap.begin(); ap != ap.end(); ap++) {
    (*ap).x = (*ap).positionBox->open();
    (*ap).f = (*ap).forceBox->open();
    (*ap).a = (*ap).atomBox->open();
  } 

  // take triplet and pass with tuple info to force eval
  DebugM(2, "ComputeHomeTuples::doWork() - size of tuple list = " << tupleList.size() << endl );
  ResizeArrayIter<T> al(tupleList);
  for (al = al.begin(); al != al.end(); al++ ) {
    // computeForce returns (BigReal)change in energy.  This must be used.
    (*al).computeForce();
  }

  for (ap = ap.begin(); ap != ap.end(); ap++) {
    (*ap).positionBox->close(&(*ap).x);
    (*ap).forceBox->close(&(*ap).f);
    (*ap).atomBox->close(&(*ap).a);
  }
  DebugM(1, "ComputeHomeTuples::doWork() -- done" << endl);
}

