
#ifndef MIGRATION_H
#define MIGRATION_H

#include "NamdTypes.h"

struct MigrationElem {
  PatchID src, dest;
  AtomID atomID;
  AtomProperties atomProp;
  Position posInit, pos;
  Velocity vel;
  Force force, forceShort, forceLong;
  int xdev, ydev, zdev;
  MigrationElem(PatchID s, AtomID aid, AtomProperties ap, Position pInit,
		Position p, Velocity v, 
		Force f, Force f_short, Force f_long, 
		int xd, int yd, int zd) : 
      src(s), atomID(aid), atomProp(ap), posInit(pInit), pos(p), 
      vel(v), force(f), forceShort(f_short), forceLong(f_long),
      xdev(xd), ydev(yd), zdev(zd) {};
}
  
typedef ResizeArray<MigrationElem> MigrationList;

#endif // MIGRATION_H
