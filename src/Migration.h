/**
***  Copyright (c) 1995, 1996, 1997, 1998, 1999, 2000 by
***  The Board of Trustees of the University of Illinois.
***  All rights reserved.
**/

#ifndef MIGRATION_H
#define MIGRATION_H

#include "InfoStream.h"
#include "NamdTypes.h"
#include "PatchTypes.h"

struct MigrationElem {
  AtomID atomID;
  AtomProperties atomProp;
  Transform trans;
  Position pos;
  Position pos_checkpoint;
  Velocity vel;
  Force force[Results::maxNumForces];
  MigrationElem() {};
  MigrationElem(AtomID &aid, AtomProperties &ap, Transform &t,
		Position &p, Position &p_c, Velocity &v, 
		Force (&f)[Results::maxNumForces]) : 
      atomID(aid), atomProp(ap), trans(t), pos(p), pos_checkpoint(p_c), vel(v)
  {
    for ( int i = 0; i < Results::maxNumForces; ++i ) force[i] = f[i];
  }
  MigrationElem(const MigrationElem &other)
  {
    atomID = other.atomID;
    atomProp = other.atomProp;
    trans = other.trans;
    pos = other.pos;
    pos_checkpoint = other.pos_checkpoint;
    vel = other.vel;
    for ( int i = 0; i < Results::maxNumForces; ++i )
      force[i] = other.force[i];
  }
  MigrationElem& operator=(const MigrationElem &other)
  {
    atomID = other.atomID;
    atomProp = other.atomProp;
    trans = other.trans;
    pos = other.pos;
    pos_checkpoint = other.pos_checkpoint;
    vel = other.vel;
    for ( int i = 0; i < Results::maxNumForces; ++i )
      force[i] = other.force[i];
    return *this;
  }
  ~MigrationElem() {}
  int operator==(const MigrationElem &m) {
    return (atomID == m.atomID);
  }
  void print() {
    iout << "Atom ID = " << atomID << endi;
    iout << "  position = " << pos << endi;
    iout << "  velocity = " << vel << endi;
    iout << "  force = " << force << endi;
  }
};

typedef ResizeArray<MigrationElem> MigrationList;
typedef ResizeArrayIter<MigrationElem> MigrationListIter;

struct MigrationInfo {
  PatchID destPatchID;
  NodeID  destNodeID;
  MigrationList mList;
};

#endif // MIGRATION_H

