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

template <class T>
ComputeHomeTuples<T>::ComputeHomeTuples(ComputeID c) : Compute(c) {
  CPrintf("ComputeHomeTuples::ComputeHomeTuples(%d) -- starting\n",(int)c);
  patchMap = PatchMap::Object();
  atomMap = AtomMap::Object();

  maxProxyAtoms = 0;
  dummy = NULL;
  CPrintf("ComputeHomeTuples::ComputeHomeTuples(%d) -- done\n",(int)c);
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
  CPrintf("ComputeHomeTuples::mapReady() - Starting Up\n");
  HomePatchList *a = patchMap->homePatchList();
  ResizeArrayIter<HomePatchElem> ai(*a);

  tuplePatchList.resize(0);
  CPrintf("ComputeHomeTuples::mapReady() - Size of the tuplePatchList %d\n",
	tuplePatchList.size());

  for ( ai = ai.begin(); ai != ai.end(); ai++ ) {
    tuplePatchList.add(TuplePatchElem((*ai).p, HOME, cid));
    CPrintf("ComputeHomeTuples::mapReady() - adding Patch %d to list\n",
      (*ai).p->getPatchID() );
  }

  /* cycle through each patch */
  CPrintf("ComputeHomeTuples::mapReady() - iterating over patches to get atoms\n");
  for ( ai = ai.begin(); ai != ai.end(); ai++ )
  {
    Patch *p = (*ai).p;
    CPrintf("ComputeHomeTuples::mapReady() - looking at patch %d\n", 
      (*ai).p->getPatchID() );
    AtomIDList &atomID = p->getAtomIDList();

    /* cycle through each angle in the patch */
    CPrintf("ComputeAtoms::mapReady() - patch has %d atoms\n", 
      p->getNumAtoms() );
    for (int i=0; i < p->getNumAtoms(); i++)
    {
      T::addTuplesForAtom((void*)&tupleList,atomID[i],node->molecule);
    }
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
  CPrintf("ComputeHomeTuples::doWork() -- started\n");
  // Open Boxes
  ResizeArrayIter<TuplePatchElem> ap(tuplePatchList);
  for (ap = ap.begin(); ap != ap.end(); ap++) {
    (*ap).x = (*ap).positionBox->open();
    (*ap).f = (*ap).forceBox->open();
  } 

  // take triplet and pass with tuple info to force eval
  ResizeArrayIter<T> al(tupleList);
  for (al = al.begin(); al != al.end(); al++ ) {
    // computeForce returns (BigReal)change in energy.  This must be used.
    (*al).computeForce();
  }

  for (ap = ap.begin(); ap != ap.end(); ap++) {
    (*ap).positionBox->close(&(*ap).x);
    (*ap).forceBox->close(&(*ap).f);
  }
  CPrintf("ComputeHomeTuples::doWork() -- done\n");
}

