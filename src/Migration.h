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
  int xdev, ydev, zdev;
  MigrationElem() {};
  MigrationElem(AtomID aid, AtomProperties ap, Position pInit,
		Position p, Velocity v, 
		Force f, Force f_short, Force f_long, 
		int xd, int yd, int zd) : 
      atomID(aid), atomProp(ap), posInit(pInit), pos(p), 
      vel(v), force(f), forceShort(f_short), forceLong(f_long),
      xdev(xd), ydev(yd), zdev(zd) {};
  int operator==(const MigrationElem &m) {
    return (atomID == m.atomID);
  }
};

typedef ResizeArray<MigrationElem> MigrationList;
typedef ResizeArrayIter<MigrationElem> MigrationListIter;

#endif // MIGRATION_H
