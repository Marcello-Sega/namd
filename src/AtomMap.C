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
 *	$Revision: 1.3 $	$Date: 1996/11/22 01:43:32 $
 *
 ***************************************************************************
 * REVISION HISTORY:
 *
 * $Log: AtomMap.C,v $
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

