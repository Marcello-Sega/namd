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

/***************************************************************************
 * RCS INFORMATION:
 *
 *	$RCSfile: AtomMap.C,v $
 *	$Author: brunner $	$Locker:  $		$State: Exp $
 *	$Revision: 1.1010 $	$Date: 1999/05/11 23:56:11 $
 *
 ***************************************************************************
 * REVISION HISTORY:
 *
 * $Log: AtomMap.C,v $
 * Revision 1.1010  1999/05/11 23:56:11  brunner
 * Changes for new charm version
 *
 * Revision 1.1009  1998/03/03 23:04:59  brunner
 * Changed include files for new simplified Charm++ include file structure.
 *
 * Revision 1.1008  1997/11/07 20:17:30  milind
 * Made NAMD to run on shared memory machines.
 *
 * Revision 1.1007  1997/04/10 09:13:46  ari
 * Final debugging for compute migration / proxy creation for load balancing.
 * Lots of debug code added, mostly turned off now.
 * Fixed bug in PositionBox when Patch had no dependencies.
 * Eliminated use of cout and misuse of iout in numerous places.
 *                                            Ari & Jim
 *
 * Revision 1.1006  1997/03/04 22:37:03  ari
 * Clean up of code.  Debug statements removal, dead code removal.
 * Minor fixes, output fixes.
 * Commented some code from the top->down.  Mainly reworked Namd, Node, main.
 *
 * Revision 1.1005  1997/02/14 19:07:27  nealk
 * Added new/delete comments.
 * Played with DPMTA.
 *
 * Revision 1.1004  1997/02/13 17:06:19  jim
 * Turned off debugging.
 *
 * Revision 1.1003  1997/02/13 16:17:09  ari
 * Intermediate debuging commit - working to fix deep bug in migration?
 *
 * Revision 1.1002  1997/02/07 16:49:26  jim
 * Fixing bugs that affect parallel atom migration.
 *
 * Revision 1.1001  1997/02/07 05:42:28  ari
 * Some bug fixing - atom migration on one node works
 * Atom migration on multiple nodes gets SIGSEGV
 *
 * Revision 1.1000  1997/02/06 15:57:27  ari
 * Resetting CVS to merge branches back into the main trunk.
 * We will stick to main trunk development as suggested by CVS manual.
 * We will set up tags to track fixed points of development/release
 * as suggested by CVS manual - all praise the CVS manual.
 *
 * Revision 1.778  1997/01/28 00:29:46  ari
 * internal release uplevel to 1.778
 *
 * Revision 1.777.2.2  1997/01/17 19:56:26  ari
 * Another Test
 *
 * Revision 1.777.2.1  1997/01/17 19:54:47  ari
 * Test for development-0-1 branch cvs commit
 *
 * Revision 1.777  1997/01/17 19:35:17  ari
 * Internal CVS leveling release.  Start development code work
 * at 1.777.1.1.
 *
 * Revision 1.3  1996/11/22 01:43:32  jim
 * switched to use AtomIDList
 *
 * Revision 1.2  1996/10/29 23:36:13  ari
 * *** empty log message ***
 *
 * Revision 1.1  1996/10/24 18:50:33  brunner
 * Initial revision
 *
 *
 ***************************************************************************/

