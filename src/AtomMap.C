/***************************************************************************/
/*                                                                         */
/*              (C) Copyright 1996 The Board of Trustees of the            */
/*                          University of Illinois                         */
/*                           All Rights Reserved                           */
/*									   */
/***************************************************************************/

/***************************************************************************
 * DESCRIPTION: Maps Atoms on node.  Singleton
 *
 ***************************************************************************/

#include "AtomMap.h"

AtomMap *AtomMap::_instance = 0;

AtomMap *AtomMap::Instance() {
  if (_instance == 0) {
    _instance = new AtomMap;
  }
  return _instance;
}

//----------------------------------------------------------------------
AtomMap::AtomMap(void)
{
  localIDTable = NULL;
  cleared = false;
}

//----------------------------------------------------------------------
AtomMap::~AtomMap(void)
{
  delete [] localIDTable;  // Delete on a NULL pointer should be ok
}

//----------------------------------------------------------------------
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
	localIDTable[al[i]].pid = pid;
	localIDTable[al[i]].index = i;
    }
    cleared = false;
    return 0;
  }
}

//----------------------------------------------------------------------
void AtomMap::clearMap(void)
{
  if (!cleared && localIDTable != NULL)
  {
    for(int i=0; i < tableSz; i++)
      localIDTable[i].pid = localIDTable[i].index = notUsed;
    cleared = true;
  }
}

/***************************************************************************
 * RCS INFORMATION:
 *
 *	$RCSfile: AtomMap.C,v $
 *	$Author: jim $	$Locker:  $		$State: Exp $
 *	$Revision: 1.1002 $	$Date: 1997/02/07 16:49:26 $
 *
 ***************************************************************************
 * REVISION HISTORY:
 *
 * $Log: AtomMap.C,v $
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

