/***************************************************************************/
/*                                                                         */
/*              (C) Copyright 1996 The Board of Trustees of the            */
/*                          University of Illinois                         */
/*                           All Rights Reserved                           */
/*									   */
/***************************************************************************/

/***************************************************************************
 * DESCRIPTION: Tracks location of Atoms on node.  Singleton.
 *
 ***************************************************************************/

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
int AtomMap::unregisterIDs(PatchID pid, AtomIDList al)
{
  if (localIDTable == NULL)
    return -1;
  else 
  {
    int ali;

    for(int i = 0; i < al.size(); ++i)
    {
	if (localIDTable[ali = al[i]].pid == pid) {
	    localIDTable[ali].pid = notUsed;
	    localIDTable[ali].index = notUsed;
	} else {
	    DebugM(4, "Avoided overwriting atomID " << ali << "\n");
	}
    }
    return 0;
  }
}
//----------------------------------------------------------------------
int AtomMap::registerIDs(PatchID pid, AtomIDList al)
{
  if (localIDTable == NULL)
    return -1;
  else 
  {
    for(int i = 0; i < al.size(); ++i)
    {
	if (localIDTable[al[i]].pid != notUsed) {
	  DebugM(4, "Overwriting atomID " << al[i] <<" used to be on " 
	    << localIDTable[al[i]].pid << " now is on " << pid << "\n");
	}
	localIDTable[al[i]].pid = pid;
	localIDTable[al[i]].index = i;
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


