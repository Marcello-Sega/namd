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

#include "common.h"
#include "structures.h"
#include "NamdTypes.h"
#include "Node.h"
#include "PatchMap.h"
#include "AtomMap.h"
#include "ComputeAngles.h"
#include "PatchMgr.h"
#include "Molecule.h"

ComputeAngles::ComputeAngles(ComputeID c) : Compute(c) {
  patchMap = PatchMap::Object();
  atomMap = AtomMap::Object();

  maxProxyAtoms = 0;
  dummy = NULL;
}

void ComputeAngles::mapAtoms() {

  // AnglePatchList contribution for proxies should be
  // gathered here.
  // Gather all home patches
  /*
  ProxyPatchList *pp = patchMap->proxyPatchList();
  ProxyPatchListIter ppi(*pp);

  maxProxyAtoms = 0;
  delete[] dummy;

  for ( ppi = ppi.begin(); ppi != ppi.end(); ppi++ ) {
    anglePatchList.add(AnglePatchElem((*ppi).patch, HOME, cid));
    if ((*ppi).patch->getNumAtoms() > maxProxyAtoms) {
      maxProxyAtoms = (*ppi).patch->getNumAtoms();
    }
  }
  dummy = new Force[maxProxyAtoms];
  */


  // Gather all home patches
  HomePatchList *a = patchMap->homePatchList();
  ResizeArrayIter<HomePatchElem> ai(*a);

  anglePatchList.resize(0);

  for ( ai = ai.begin(); ai != ai.end(); ai++ ) {
    anglePatchList.add(AnglePatchElem((*ai).p, HOME, cid));
  }

  /* cycle through each patch */
  for ( ai = ai.begin(); ai != ai.end(); ai++ )
  {
    Patch *p = (*ai).p;
    AtomIDList &atomID = p->getAtomIDList();

    /* cycle through each angle in the patch */
    for (int i=0; i < p->getNumAtoms(); i++)
    {
      /* get list of all angles for the atom */
      LintList *angles = node->molecule->get_angles_for_atom(atomID[i]);

      /* cycle through each angle */
      int angleNum = angles->head();
      while(angleNum != LIST_EMPTY)
      {
        /* store angle in the list */
        angleList.add(AngleElem(node->molecule->get_angle(angleNum)));
        angleNum = angles->next();
      }
    }
  }

  // Resolve all atoms in angleList to correct PatchList element and index
  ResizeArrayIter<AngleElem> al(angleList);

  for (al = al.begin(); al != al.end(); al++ ) {
    for (int i=0; i < 3; i++) {
	LocalID aid = atomMap->localID((*al).atomID[i]);
	(*al).p[i] = anglePatchList.find(AnglePatchElem(aid.pid));
	(*al).localIndex[i] = aid.index;
    }
  }
}


void ComputeAngles::doWork() {
  // Open Boxes
  ResizeArrayIter<AnglePatchElem> ap(anglePatchList);
  for (ap = ap.begin(); ap != ap.end(); ap++) {
    (*ap).x = (*ap).positionBox->open();
    (*ap).f = (*ap).forceBox->open();
  } 

  // take triplet and pass with angle info to force eval
  ResizeArrayIter<AngleElem> al(angleList);
  for (al = al.begin(); al != al.end(); al++ ) {
    angleForce((*al).p[0]->x[(*al).localIndex[0]],
	       (*al).p[1]->x[(*al).localIndex[1]],
	       (*al).p[2]->x[(*al).localIndex[2]],
	       (*al).p[0]->f+(*al).localIndex[0],
	       (*al).p[1]->f+(*al).localIndex[1],
	       (*al).p[2]->f+(*al).localIndex[2],
	       (*al).angleType);
  }

  for (ap = ap.begin(); ap != ap.end(); ap++) {
    (*ap).positionBox->close(&(*ap).x);
    (*ap).forceBox->close(&(*ap).f);
  }
}
