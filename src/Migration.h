//-*-c++-*-

#ifndef MIGRATION_H
#define MIGRATION_H

#include "InfoStream.h"
#include "NamdTypes.h"
#include "PatchTypes.h"

struct MigrationElem {
  AtomID atomID;
  AtomProperties atomProp;
  Position posInit, pos;
  Velocity vel;
  Force force[Results::maxNumForces];
  MigrationElem() {};
  MigrationElem(AtomID aid, AtomProperties ap, Position pInit,
		Position p, Velocity v, 
		Force f[Results::maxNumForces]) : 
      atomID(aid), atomProp(ap), posInit(pInit), pos(p), vel(v)
  {
    for ( int i = 0; i < Results::maxNumForces; ++i ) force[i] = f[i];
  }
  MigrationElem(const MigrationElem &other)
  {
    atomID = other.atomID;
    atomProp = other.atomProp;
    posInit = other.posInit;
    pos = other.pos;
    vel = other.vel;
    for ( int i = 0; i < Results::maxNumForces; ++i )
      force[i] = other.force[i];
  }
  MigrationElem& operator=(const MigrationElem &other)
  {
    atomID = other.atomID;
    atomProp = other.atomProp;
    posInit = other.posInit;
    pos = other.pos;
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
  MigrationList *mList;
};

#endif // MIGRATION_H
