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

#ifndef ATOMMAP_H
#define ATOMMAP_H

#include "NamdTypes.h"

enum { notUsed = -1 };

class AtomMap
{
public:
  static AtomMap *Instance();
  inline static AtomMap *Object() {return CpvAccess(AtomMap_instance);}
  ~AtomMap(void);
  void checkMap();

  void allocateMap(int nAtomIDs);

  int registerIDs(PatchID pid, AtomIDList al);
  int unregisterIDs(PatchID pid, AtomIDList al);

  LocalID localID(AtomID id);

  void clearMap(void);
  void print(void);


protected:
  AtomMap(void);

private:

  LocalID *localIDTable;
  int tableSz;
  Boolean cleared;

};

//----------------------------------------------------------------------
// LocalID contains patch pid and local patch atom index
// for a given global atom number
inline LocalID AtomMap::localID(AtomID id)
{
  return localIDTable[id];
}

#endif /* ATOMMAP_H */

