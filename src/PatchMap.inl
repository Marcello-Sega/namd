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
inline int PatchMap::numPatches(void)
{
  return nPatches;
}

//----------------------------------------------------------------------
inline PatchID PatchMap::assignToPatch(Position p)
{
  int xi, yi, zi;
  xi = (int)floor(((BigReal)xDim)*((p.x-xOrigin)/xLength));
  yi = (int)floor(((BigReal)yDim)*((p.y-yOrigin)/yLength));
  zi = (int)floor(((BigReal)zDim)*((p.z-zOrigin)/zLength));
  return pid(xi,yi,zi);
}

//----------------------------------------------------------------------
inline int PatchMap::xDimension(void)
{
  return xDim;
}

//----------------------------------------------------------------------
inline int PatchMap::yDimension(void)
{
  return yDim;
}

//----------------------------------------------------------------------
inline int PatchMap::zDimension(void)
{
  return zDim;
}

//----------------------------------------------------------------------
inline Vector PatchMap::Origin(void)
{
  Vector o;
  o.x = xOrigin;
  o.y = yOrigin;
  o.z = zOrigin;
  return(o);
}

//----------------------------------------------------------------------
inline int PatchMap::xIsPeriodic(void)
{
  return xPeriodic;
}

//----------------------------------------------------------------------
inline int PatchMap::yIsPeriodic(void)
{
  return yPeriodic;
}

//----------------------------------------------------------------------
inline int PatchMap::zIsPeriodic(void)
{
  return zPeriodic;
}

//----------------------------------------------------------------------
#define MODULO(I,J) ( (I)<0 ? (I)-(J)*((I)/(J)-1) : (I)-(J)*((I)/(J)) )

inline int PatchMap::pid(int xIndex, int yIndex, int zIndex)
{
  int allsame = 0;
  if ( xPeriodic ) xIndex = MODULO(xIndex,xDim);
  else
  {
    if ( xIndex < 0 ) xIndex = 0;
    if ( xIndex >= xDim ) xIndex = xDim - 1;
  }
  if ( yPeriodic ) yIndex = MODULO(yIndex,yDim);
  else
  {
    if ( yIndex < 0 ) yIndex = 0;
    if ( yIndex >= yDim ) yIndex = yDim - 1;
  }
  if ( zPeriodic ) zIndex = MODULO(zIndex,zDim);
  else
  {
    if ( zIndex < 0 ) zIndex = 0;
    if ( zIndex >= zDim ) zIndex = zDim - 1;
  }
  return ((zIndex*yDim)+yIndex)*xDim + xIndex;
}

//----------------------------------------------------------------------
inline int PatchMap::xIndex(int pid)
{
  return pid % xDim;
}

//----------------------------------------------------------------------
inline int PatchMap::yIndex(int pid)
{
  return (pid / xDim) % yDim;
}

//----------------------------------------------------------------------
inline int PatchMap::zIndex(int pid)
{
  return (pid / (xDim*yDim));
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

    // z
    register int k = pdat1->zi;
    register int k2 = pdat2->zi;
    if ( ( k ? k : zMaxIndex ) == k2 + 1 ) k = k2;

    // y
    register int j = pdat1->yi;
    register int j2 = pdat2->yi;
    if ( ( j ? j : yMaxIndex ) == j2 + 1 ) j = j2;

    // x
    register int i = pdat1->xi;
    register int i2 = pdat2->xi;
    if ( ( i ? i : xMaxIndex ) == i2 + 1 ) i = i2;

    ds = ((k*yDim)+j)*xDim + i;
  }

  return ds;
}

//----------------------------------------------------------------------
inline int PatchMap::node(int pid)
{
  return patchData[pid].node;
}

//----------------------------------------------------------------------
inline Coordinate PatchMap::minX(int pid)
{
  return patchData[pid].x0;
}

//----------------------------------------------------------------------
inline Coordinate PatchMap::maxX(int pid)
{
  return patchData[pid].x1;
}

//----------------------------------------------------------------------
inline Coordinate PatchMap::minY(int pid)
{
  return patchData[pid].y0;
}

//----------------------------------------------------------------------
inline Coordinate PatchMap::maxY(int pid)
{
  return patchData[pid].y1;
}

//----------------------------------------------------------------------
inline Coordinate PatchMap::minZ(int pid)
{
  return patchData[pid].z0;
}

//----------------------------------------------------------------------
inline Coordinate PatchMap::maxZ(int pid)
{
  return patchData[pid].z1;
}

//----------------------------------------------------------------------
inline int PatchMap::numCids(int pid)
{
  return patchData[pid].numCids;
}

//----------------------------------------------------------------------
inline int PatchMap::cid(int pid,int i)
{
  return patchData[pid].cids[i];
}

#endif
/***************************************************************************
 * RCS INFORMATION:
 *
 *	$RCSfile: PatchMap.inl,v $
 *	$Author: jim $	$Locker:  $		$State: Exp $
 *	$Revision: 1.2 $	$Date: 1998/07/16 18:52:14 $
 *
 ***************************************************************************
 * REVISION HISTORY:
 *
 * $Log: PatchMap.inl,v $
 * Revision 1.2  1998/07/16 18:52:14  jim
 * Localized common downstream patch optimization.
 *
 * Revision 1.1  1997/10/06 00:12:34  jim
 * Added PatchMap.inl, sped up cycle-boundary tuple code.
 *
 *
 ***************************************************************************/

