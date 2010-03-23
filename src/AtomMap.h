/**
***  Copyright (c) 1995, 1996, 1997, 1998, 1999, 2000 by
***  The Board of Trustees of the University of Illinois.
***  All rights reserved.
**/

#ifndef ATOMMAP_H
#define ATOMMAP_H

#include "NamdTypes.h"

enum { notUsed = -1 };

class AtomMap
{
public:
  static AtomMap *Instance();
  inline static AtomMap *Object() {return CkpvAccess(AtomMap_instance);}
  ~AtomMap(void);
  void checkMap();

  void allocateMap(int nAtomIDs);

  int registerIDs(PatchID pid, const CompAtomExt *begin, const CompAtomExt *end);
  int registerIDsFullAtom(PatchID pid, const FullAtom *begin, const FullAtom *end);
  int unregisterIDs(PatchID pid, const CompAtomExt *begin, const CompAtomExt *end);
  int unregisterIDs(PatchID pid, const FullAtom *begin, const FullAtom *end);

  LocalID localID(AtomID id);

  void clearMap(void);
  void print(void);


protected:
  AtomMap(void);

private:

  LocalID *localIDTable;
  int tableSz;
  bool cleared;

};

//----------------------------------------------------------------------
// LocalID contains patch pid and local patch atom index
// for a given global atom number
inline LocalID AtomMap::localID(AtomID id)
{
  return localIDTable[id];
}

#endif /* ATOMMAP_H */

