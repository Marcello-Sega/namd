//-*-c++-*-

#ifndef MIGRATION_H
#define MIGRATION_H

#include "NamdTypes.h"


struct MigrationElem {
  AtomID atomID;
  AtomProperties atomProp;
  Position posInit, pos;
  Velocity vel;
  Force force, forceShort, forceLong;
  MigrationElem() {};
  MigrationElem(AtomID aid, AtomProperties ap, Position pInit,
		Position p, Velocity v, 
		Force f, Force f_short, Force f_long) : 
      atomID(aid), atomProp(ap), posInit(pInit), pos(p), 
      vel(v), force(f), forceShort(f_short), forceLong(f_long) {};
  ~MigrationElem() {}
  int operator==(const MigrationElem &m) {
    return (atomID == m.atomID);
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
