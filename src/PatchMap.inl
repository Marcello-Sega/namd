/***************************************************************************/
/*        (C) Copyright 1996,1997 The Board of Trustees of the             */
/*                          University of Illinois                         */
/*                           All Rights Reserved                           */
/***************************************************************************/
/***************************************************************************
 * DESCRIPTION:
 *  Inline functions for PatchMap.
 *
 ***************************************************************************/

#ifndef PATCHMAP_INL
#define PATCHMAP_INL

#include "PatchMap.h"
#include "AtomMap.h"

//----------------------------------------------------------------------
inline PatchID PatchMap::assignToPatch(Position p, const Lattice &l)
{
  int ai, bi, ci;
  ScaledPosition s = l.scale(p);
  ai = (int)floor(((BigReal)aDim)*((s.x-aOrigin)/aLength));
  bi = (int)floor(((BigReal)bDim)*((s.y-bOrigin)/bLength));
  ci = (int)floor(((BigReal)cDim)*((s.z-cOrigin)/cLength));
  return pid(ai,bi,ci);
}

//----------------------------------------------------------------------
#define MODULO(I,J) ( (I)<0 ? ((J)-(-1*(I))%(J))%(J) : (I)%(J) )

inline int PatchMap::pid(int aIndex, int bIndex, int cIndex)
{
  if ( aPeriodic ) aIndex = MODULO(aIndex,aDim);
  else
  {
    if ( aIndex < 0 ) aIndex = 0;
    if ( aIndex >= aDim ) aIndex = aDim - 1;
  }
  if ( bPeriodic ) bIndex = MODULO(bIndex,bDim);
  else
  {
    if ( bIndex < 0 ) bIndex = 0;
    if ( bIndex >= bDim ) bIndex = bDim - 1;
  }
  if ( cPeriodic ) cIndex = MODULO(cIndex,cDim);
  else
  {
    if ( cIndex < 0 ) cIndex = 0;
    if ( cIndex >= cDim ) cIndex = cDim - 1;
  }
  return ((cIndex*bDim)+bIndex)*aDim + aIndex;
}

//----------------------------------------------------------------------
inline int PatchMap::downstream(int pid1, int pid2)
{
  register int ds;

  if ( pid1 == pid2 ) { ds = pid1; }

  else if ( pid1 == notUsed || pid2 == notUsed ) { ds =  notUsed; }

  else {
    register PatchData *pdat1 = &(patchData[pid1]);
    register PatchData *pdat2 = &(patchData[pid2]);

    // c
    register int k = pdat1->cIndex;
    register int k2 = pdat2->cIndex;
    if ( ( k ? k : cMaxIndex ) == k2 + 1 ) k = k2;

    // b
    register int j = pdat1->bIndex;
    register int j2 = pdat2->bIndex;
    if ( ( j ? j : bMaxIndex ) == j2 + 1 ) j = j2;

    // a
    register int i = pdat1->aIndex;
    register int i2 = pdat2->aIndex;
    if ( ( i ? i : aMaxIndex ) == i2 + 1 ) i = i2;

    ds = ((k*bDim)+j)*aDim + i;
  }

  return ds;
}

#endif
/***************************************************************************
 * RCS INFORMATION:
 *
 *	$RCSfile: PatchMap.inl,v $
 *	$Author: jim $	$Locker:  $		$State: Exp $
 *	$Revision: 1.4 $	$Date: 1999/09/03 20:46:19 $
 *
 ***************************************************************************
 * REVISION HISTORY:
 *
 * $Log: PatchMap.inl,v $
 * Revision 1.4  1999/09/03 20:46:19  jim
 * Support for non-orthogonal periodic boundary conditions.
 *
 * Revision 1.3  1998/08/17 23:29:53  jim
 * Fixed MODULO macro needed for negative arguments.  I can't do math.
 *
 * Revision 1.2  1998/07/16 18:52:14  jim
 * Localized common downstream patch optimization.
 *
 * Revision 1.1  1997/10/06 00:12:34  jim
 * Added PatchMap.inl, sped up cycle-boundary tuple code.
 *
 *
 ***************************************************************************/

