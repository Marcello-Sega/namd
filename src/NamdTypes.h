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

#ifndef NAMDTYPES_H

#define NAMDTYPES_H

#include "Vector.h"
#include "Templates/ResizeArray.h"
#include "Templates/ResizeArrayIter.h"

class Patch;
class Compute;

typedef Vector Position;
typedef Vector Velocity;
typedef Vector Force;
typedef int AtomID;

typedef double Coordinate;

typedef ResizeArray<Position> PositionList;
typedef ResizeArrayIter<Position> PositionListIter;
typedef ResizeArray<Velocity> VelocityList;
typedef ResizeArrayIter<Velocity> VelocityListIter;
typedef ResizeArray<Force> ForceList;
typedef ResizeArrayIter<Force> ForceListIter;

typedef ResizeArray<AtomID> AtomIDList;
typedef ResizeArrayIter<AtomID> AtomIDListIter;

typedef int PatchID;
typedef int ComputeID;
typedef int NodeID;

typedef ResizeArray<PatchID> PatchIDList;
typedef ResizeArrayIter<PatchID> PatchIDListIter;
typedef ResizeArray<Patch *> PatchList;
typedef ResizeArrayIter<Patch *> PatchListIter;

typedef UniqueSortedArray<ComputeID> ComputeIDList;
typedef ResizeArrayIter<ComputeID> ComputeIDListIter;

typedef ResizeArray<Compute *> ComputeList;
typedef ResizeArrayIter<Compute *> ComputeListIter;

enum ComputeType
{
  electForceType,
  bondForceType,
  angleForceType,
  dihedralForceType,
  improperForceType
};

enum Boolean
{
  false=0,
  true=1
};

#endif /* NAMDTYPES_H */

/***************************************************************************
 * RCS INFORMATION:
 *
 *	$RCSfile: NamdTypes.h,v $
 *	$Author: ari $	$Locker:  $		$State: Exp $
 *	$Revision: 1.8 $	$Date: 1996/10/16 08:22:39 $
 *
 ***************************************************************************
 * REVISION HISTORY:
 *
 * $Log: NamdTypes.h,v $
 * Revision 1.8  1996/10/16 08:22:39  ari
 * *** empty log message ***
 *
 * Revision 1.7  1996/09/10 04:16:28  jim
 * Added iterators for all lists.
 *
 * Revision 1.6  1996/08/29 00:52:06  ari
 * *** empty log message ***
 *
 * Revision 1.5  1996/08/23 21:36:58  ari
 * *** empty log message ***
 *
 * Revision 1.4  1996/08/19 21:37:02  brunner
 * Added Coordinate
 *
 * Revision 1.3  1996/08/19 21:27:51  ari
 * .
 *
 * Revision 1.2  1996/08/16 21:42:58  brunner
 * *** empty log message ***
 *
 * Revision 1.1  1996/08/16 21:00:37  brunner
 * Initial revision
 *
 ***************************************************************************/
