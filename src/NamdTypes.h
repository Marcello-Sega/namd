/**
***  Copyright (c) 1995, 1996, 1997, 1998, 1999, 2000 by
***  The Board of Trustees of the University of Illinois.
***  All rights reserved.
**/

#ifndef NAMDTYPES_H

#define NAMDTYPES_H

#include "Vector.h"
#include "ResizeArray.h"

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
  char i,j,k;
  Transform(void) { i=0; j=0; k=0; }
};

struct CompAtom {
  Position position;
  Charge charge;
  unsigned int id : 24;
  unsigned int hydrogenGroupSize : 3;
  unsigned int nonbondedGroupSize : 3;
  unsigned int atomFixed : 1;
  unsigned int groupFixed : 1;

  CompAtom() { ; }

  // Needed for IBM's xlC compiler
  inline CompAtom(const CompAtom &a) :
    position(a.position), charge(a.charge),
    id(a.id), hydrogenGroupSize(a.hydrogenGroupSize),
    nonbondedGroupSize(a.nonbondedGroupSize),
    atomFixed(a.atomFixed), groupFixed(a.groupFixed) {
    ;
  }

  // Needed for IBM's xlC compiler
  inline CompAtom& operator=(const CompAtom &a) {
    position = a.position;
    charge = a.charge;
    id = a.id;
    hydrogenGroupSize = a.hydrogenGroupSize;
    nonbondedGroupSize = a.nonbondedGroupSize;
    atomFixed = a.atomFixed;
    groupFixed = a.groupFixed;
    return *this;
  }

};

struct FullAtom : CompAtom {
  Velocity velocity;
  Mass mass;
  Transform transform;
};

typedef ResizeArray<CompAtom> CompAtomList;
typedef ResizeArray<FullAtom> FullAtomList;
typedef ResizeArray<Position> PositionList;
typedef ResizeArray<Velocity> VelocityList;
typedef ResizeArray<Force> ForceList;
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

