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
#include "ComputeHomePatches.h"
#include "PatchMgr.h"
#include "HomePatchList.h"
#include "Molecule.h"
#include "ReductionMgr.h"
#define MIN_DEBUG_LEVEL 3
#include "Debug.h"

ComputeHomePatches::ComputeHomePatches(ComputeID c) : Compute(c) {
  patchMap = PatchMap::Object();
  reduction = ReductionMgr::Object();

  fake_seq = 0;
}

ComputeHomePatches::~ComputeHomePatches()
{
  ;
}

void ComputeHomePatches::initialize()
{
  HomePatchList *a = patchMap->homePatchList();
  ResizeArrayIter<HomePatchElem> ai(*a);

  patchList.resize(0);

  for ( ai = ai.begin(); ai != ai.end(); ai++ ) {
    patchList.add(PatchElem((*ai).patch, cid));
  }

  setNumPatches(patchList.size());
}

void ComputeHomePatches::atomUpdate()
{
  Compute::atomUpdate();
}

