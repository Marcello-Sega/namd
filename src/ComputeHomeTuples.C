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

#define  MIN_DEBUG_LEVEL 3
#define  DEBUGM
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

  // ComputeHomeTuples contribution for proxies should be
  // gathered here.
  // Gather all home patches
  /*
  ProxyPatchList *pp = patchMap->proxyPatchList();
  ProxyPatchListIter ppi(*pp);

  maxProxyAtoms = 0;
  delete[] dummy;

  for ( ppi = ppi.begin(); ppi != ppi.end(); ppi++ ) {
    tuplePatchList.add(TuplePatchElem((*ppi).patch, HOME, cid));
    if ((*ppi).patch->getNumAtoms() > maxProxyAtoms) {
      maxProxyAtoms = (*ppi).patch->getNumAtoms();
    }
  }
  dummy = new Force[maxProxyAtoms];
  */


  // Gather all home patches
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

  /* cycle through each patch */
  DebugM(1, "ComputeHomeTuples::mapReady() - iterating over patches to get atoms" << endl);
  for ( ai = ai.begin(); ai != ai.end(); ai++ )
  {
    Patch *p = (*ai).p;
    DebugM(1, "ComputeHomeTuples::mapReady() - looking at patch " << (*ai).p->getPatchID() << endl );
    AtomIDList &atomID = p->getAtomIDList();

    /* cycle through each angle in the patch */
    DebugM(1, "ComputeHomeTuples::mapReady() - patch has " << p->getNumAtoms() << " atoms\n" );

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
  } 

  // take triplet and pass with tuple info to force eval
  DebugM(1, "ComputeHomeTuples::doWork() - size of atom list = " << tupleList.size() << endl );
  ResizeArrayIter<T> al(tupleList);
  for (al = al.begin(); al != al.end(); al++ ) {
    // computeForce returns (BigReal)change in energy.  This must be used.
    DebugM(4,"Atom 0 is " << (*al).p[0]->x[0] << "\n");
    (*al).computeForce();
  }

  for (ap = ap.begin(); ap != ap.end(); ap++) {
    (*ap).positionBox->close(&(*ap).x);
    (*ap).forceBox->close(&(*ap).f);
  }
  DebugM(1, "ComputeHomeTuples::doWork() -- done" << endl);
}

