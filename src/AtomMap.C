/**
***  Copyright (c) 1995, 1996, 1997, 1998, 1999, 2000 by
***  The Board of Trustees of the University of Illinois.
***  All rights reserved.
**/

/*
   Tracks location of Atoms on node.  Singleton.
*/

#include "charm++.h"

#include "ProcessorPrivate.h"

#include "AtomMap.h"

#define MIN_DEBUG_LEVEL 4
// #define DEBUGM
#include "Debug.h"


// Singleton method
AtomMap *AtomMap::Instance() {
  if (CpvAccess(AtomMap_instance) == 0) {
    CpvAccess(AtomMap_instance) = new AtomMap;	// this is never deleted!
  }
  return CpvAccess(AtomMap_instance);
}

//----------------------------------------------------------------------
AtomMap::AtomMap(void)
{
  localIDTable = NULL;
  cleared = false;
}

void
AtomMap::checkMap(void)
{ }
  

//----------------------------------------------------------------------
AtomMap::~AtomMap(void)
{
  delete [] localIDTable;  // Delete on a NULL pointer should be ok
}

//----------------------------------------------------------------------
// Creates fixed size table
void AtomMap::allocateMap(int nAtomIds)
{
  localIDTable = new LocalID[nAtomIds];
  tableSz = nAtomIds;
  for(int i=0; i < nAtomIds; i++)
    localIDTable[i].pid = localIDTable[i].index = notUsed;
  cleared = true;
}

//
int AtomMap::unregisterIDs(PatchID pid, const CompAtom *begin, const CompAtom *end)
{
  if (localIDTable == NULL)
    return -1;
  else 
  {
    for(const CompAtom *a = begin; a != end; ++a)
    {
        unsigned int ali = a->id;
	if (localIDTable[ali].pid == pid) {
	    localIDTable[ali].pid = notUsed;
	    localIDTable[ali].index = notUsed;
	}
    }
    return 0;
  }
}
//----------------------------------------------------------------------
int AtomMap::registerIDs(PatchID pid, const CompAtom *begin, const CompAtom *end)
{
  if (localIDTable == NULL)
    return -1;
  else 
  {
    for(const CompAtom *a = begin; a != end; ++a)
    {
        unsigned int ali = a->id;
	localIDTable[ali].pid = pid;
	localIDTable[ali].index = a - begin;
    }
    cleared = false;
    return 0;
  }
}

//----------------------------------------------------------------------
// resets map to notUsed condition
void AtomMap::clearMap(void)
{
  if (!cleared && localIDTable != NULL)
  {
    for(int i=0; i < tableSz; i++)
      localIDTable[i].pid = localIDTable[i].index = notUsed;
    cleared = true;
  }
}

void AtomMap::print()
{
  for (int i=0; i<tableSz; i++) {
    CkPrintf("AtomMap on node %d\n", CkMyPe());
    CkPrintf("AtomID %d -> PatchID %d:Index %d\n", i, localIDTable[i].pid,
      localIDTable[i].index);
  }
}


