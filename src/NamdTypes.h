/**
***  Copyright (c) 1995, 1996, 1997, 1998, 1999, 2000 by
***  The Board of Trustees of the University of Illinois.
***  All rights reserved.
**/

#ifndef NAMDTYPES_H

#define NAMDTYPES_H

#include "Vector.h"
#include "ResizeArray.h"
#include "ResizeArrayIter.h"

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

struct Transform
{
  int i,j,k;
  Transform(void) { i=0; j=0; k=0; }
  int operator==(const Transform &o) const {
	return ( i==o.i && j==o.j && k==o.k ); }
};

struct AtomProperties
{
  AtomID id;
  AtomType type;
  Mass mass;
  Charge charge;

  // other data information
  char hydrogenGroupSize;	// 0 from group members, !0 for group parents
  char nonbondedGroupSize;	// same, but variable with strict size limits
  // Bool water;	// TRUE if water atom (O or H)  NEVER USED -JCP
  unsigned char flags;	// for fixed atoms, etc. - use with & operator

  int operator==(const AtomProperties& a) {
    return( id == a.id );
  }
};

// Definitions for AtomProperties flags
#define ATOM_FIXED	0x0001
#define GROUP_FIXED	0x0002

typedef ResizeArray<Position> PositionList;
typedef ResizeArrayIter<Position> PositionListIter;
typedef ResizeArray<Velocity> VelocityList;
typedef ResizeArrayIter<Velocity> VelocityListIter;
typedef ResizeArray<Force> ForceList;
typedef ResizeArrayIter<Force> ForceListIter;
typedef ResizeArray<AtomProperties> AtomPropertiesList;
typedef ResizeArrayIter<AtomProperties> AtomPropertiesListIter;
typedef ResizeArray<Transform> TransformList;

typedef ResizeArray<AtomID> AtomIDList;

typedef int PatchID;
typedef int ComputeID;
typedef int NodeID;

typedef ResizeArray<PatchID> PatchIDList;
typedef ResizeArray<Patch *> PatchList;

typedef ResizeArray<Compute *> ComputeList;

// See AtomMap
struct LocalID
{
  PatchID pid;
  int index;
};

// HP compiler complains that true, false "Will be" future reserved words.
//enum Boolean
//{
//  false=0,
//  true=1
//};
#ifndef BOOLTYPE
typedef int Boolean;
#define false 0
#define true 1
#endif

#endif /* NAMDTYPES_H */

