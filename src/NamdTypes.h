//-*-c++-*-
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
#include "Templates/UniqueSortedArray.h"
#include "Templates/ResizeArrayIter.h"
#include "Templates/ResizeArrayPrimIter.h"

class Patch;
class Compute;

typedef Vector Position;
typedef Vector Velocity;
typedef Vector Force;
typedef int AtomID;
typedef int AtomType;
typedef float Mass;
typedef float Charge;

typedef double Coordinate;

struct AtomProperties
{
  AtomID id;
  AtomType type;
  Mass mass;
  Charge charge;

  int operator==(const AtomProperties& a) {
    return( id == a.id );
  }
};

typedef ResizeArray<Position> PositionList;
typedef ResizeArrayIter<Position> PositionListIter;
typedef ResizeArray<Velocity> VelocityList;
typedef ResizeArrayIter<Velocity> VelocityListIter;
typedef ResizeArray<Force> ForceList;
typedef ResizeArrayIter<Force> ForceListIter;
typedef ResizeArray<AtomProperties> AtomPropertiesList;
typedef ResizeArrayIter<AtomProperties> AtomPropertiesListIter;

typedef ResizeArray<AtomID> AtomIDList;

typedef int PatchID;
typedef int ComputeID;
typedef int NodeID;

typedef ResizeArray<PatchID> PatchIDList;
typedef ResizeArray<Patch *> PatchList;

typedef UniqueSortedArray<ComputeID> ComputeIDList;
typedef ResizeArrayPrimIter<ComputeID> ComputeIDListIter;

typedef ResizeArray<Compute *> ComputeList;

struct LocalID
{
  PatchID pid;
  int index;
};

enum ComputeType
{
  computeNonbondedSelfType,
  computeNonbondedPairType,
  computeNonbondedExclType,
  computeBondsType,
  computeAnglesType,
  computeDihedralsType,
  computeImpropersType,
  computeDPMTAType,
  computeDPMEType,
  computeFullDirectType
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
 *	$Revision: 1.779 $	$Date: 1997/02/06 15:53:17 $
 *
 ***************************************************************************
 * REVISION HISTORY:
 *
 * $Log: NamdTypes.h,v $
 * Revision 1.779  1997/02/06 15:53:17  ari
 * Updating Revision Line, getting rid of branches
 *
 * Revision 1.778.2.1  1997/02/06 05:00:44  jim
 * Added creation of full electrostatics objects.
 *
 * Revision 1.778  1997/01/28 00:30:58  ari
 * internal release uplevel to 1.778
 *
 * Revision 1.777.2.1  1997/01/27 22:45:28  ari
 * Basic Atom Migration Code added.
 * Added correct magic first line to .h files for xemacs to go to C++ mode.
 * Compiles and runs without migration turned on.
 *
 * Revision 1.777  1997/01/17 19:36:33  ari
 * Internal CVS leveling release.  Start development code work
 * at 1.777.1.1.
 *
 * Revision 1.14  1997/01/15 17:09:43  ari
 * minor changes
 *
 * Revision 1.13  1996/12/05 01:47:07  ari
 * added == to AtomProperties definition
 *
 * Revision 1.12  1996/11/30 01:27:34  jim
 * switched to realistic ComputeType definitions
 *
 * Revision 1.11  1996/10/30 01:05:50  jim
 * added AtomPropertiesList
 *
 * Revision 1.10  1996/10/30 00:32:18  jim
 * added AtomProperties definition
 *
 * Revision 1.9  1996/10/24 18:51:09  brunner
 * Added LocalID
 *
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
