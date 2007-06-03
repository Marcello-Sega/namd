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

//#ifdef ARCH_POWERPC
//typedef AlignVector Force;
//#else
typedef Vector Force;
//#endif

typedef int AtomID;
typedef int AtomType;
typedef float Mass;
typedef float Charge;

#ifdef MEM_OPT_VERSION
typedef unsigned short AtomSigID;
typedef unsigned short ExclSigID;
typedef unsigned short VDW_TYPE;
#endif

typedef double Coordinate;

struct Transform
{
  signed char i,j,k;
  Transform(void) { i=0; j=0; k=0; }
};

/*
 * 1. "position" field in this structure is very important since it
 * needs to be sent to every patch after every timestep.
 * 2. Anything that is static (value is decided before computation)
 * or only changes after atom migration should be put into the CompAtomExt structure
 * 3. This data structure is 32-byte long which is particularly optimized for some machines
 * (including BG/L) for better cache and message performance. Therefore, changes
 * to this structure should be cautioned for the sake of performance.
 */
struct CompAtom {
  Position position;
  Charge charge;
  unsigned int id : 22;
  unsigned int hydrogenGroupSize : 3;
  unsigned int nonbondedGroupIsAtom : 1;
  unsigned int atomFixed : 1;
  unsigned int groupFixed : 1;
  unsigned int partition : 4;

  CompAtom() { ; }

  // Needed for IBM's xlC compiler
  inline CompAtom(const CompAtom &a) :
    position(a.position), charge(a.charge),
    id(a.id), hydrogenGroupSize(a.hydrogenGroupSize),
    nonbondedGroupIsAtom(a.nonbondedGroupIsAtom),
    atomFixed(a.atomFixed), groupFixed(a.groupFixed),
    partition(a.partition){
      ;
  }

  // Needed for IBM's xlC compiler
  inline CompAtom& operator=(const CompAtom &a) {
    position = a.position;
    charge = a.charge;
    id = a.id;
    hydrogenGroupSize = a.hydrogenGroupSize;
    nonbondedGroupIsAtom = a.nonbondedGroupIsAtom;
    atomFixed = a.atomFixed;
    groupFixed = a.groupFixed;
    partition = a.partition;

    return *this;
  }

};

#ifdef MEM_OPT_VERSION
struct CompAtomExt {
  AtomSigID sigId;
  ExclSigID exclId;
  VDW_TYPE vdwType;

  CompAtomExt(){;}

  // Needed for IBM's xlC compiler
  inline CompAtomExt(const CompAtomExt &a) :
    sigId(a.sigId), exclId(a.exclId), vdwType(a.vdwType){
  }

  // Needed for IBM's xlC compiler
  inline CompAtomExt& operator=(const CompAtomExt &a) {
    sigId = a.sigId;
    exclId = a.exclId;
    vdwType = a.vdwType;

    return *this;
  }

};
#endif

#ifdef MEM_OPT_VERSION
struct FullAtom : CompAtom, CompAtomExt{
#else
struct FullAtom : CompAtom {
#endif
  Velocity velocity;
  Position fixedPosition;
  Mass mass;
  Transform transform;
};

typedef ResizeArray<CompAtom> CompAtomList;

#ifdef MEM_OPT_VERSION
typedef ResizeArray<CompAtomExt> CompAtomExtList;
#endif

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

typedef ResizeArray<NodeID> NodeIDList;

struct ExtForce {
  int replace;
  Force force;
  ExtForce() : replace(0) {;}
};

#endif /* NAMDTYPES_H */

