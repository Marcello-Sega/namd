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


/***************************************************************************
 * RCS INFORMATION:
 *
 *	$RCSfile $
 *	$Author $	$Locker:  $		$State: Exp $
 *	$Revision: 1.1006 $	$Date: 1999/06/17 15:46:07 $
 *
 ***************************************************************************
 * REVISION HISTORY:
 *
 * $Log: ComputeHomePatches.C,v $
 * Revision 1.1006  1999/06/17 15:46:07  jim
 * Completely rewrote reduction system to eliminate need for sequence numbers.
 *
 * Revision 1.1005  1997/09/28 10:19:04  milind
 * Fixed priorities, ReductionMgr etc.
 *
 * Revision 1.1004  1997/08/26 16:26:12  jim
 * Revamped prioritites for petter performance and easier changes.
 *
 * Revision 1.1003  1997/03/19 11:54:08  ari
 * Add Broadcast mechanism.
 * Fixed RCS Log entries on files that did not have Log entries.
 * Added some register variables to Molecule and ComputeNonbondedExcl.C
 *
 ***************************************************************************/
