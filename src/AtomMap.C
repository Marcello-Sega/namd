/***************************************************************************/
/*                                                                         */
/*              (C) Copyright 1996 The Board of Trustees of the            */
/*                          University of Illinois                         */
/*                           All Rights Reserved                           */
/*									   */
/***************************************************************************/

#include "AtomMap.h"

//----------------------------------------------------------------------
AtomMap::AtomMap(void)
{
  handle = this;
  localIDTable = NULL;
  cleared = false;
}

//----------------------------------------------------------------------
AtomMap::~AtomMap(void)
{
  handle = NULL;
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
int AtomMap::registerIDs(PatchID pid, AtomID *atomId, int nIds)
{
  if (localIDTable == NULL)
    return -1;
  else 
  {
    for(int i=0; i < nIds; i++)
    {
      if ((atomId[i] >= 0) && (atomId[i] < tableSz))
      {
	localIDTable[atomId[i]].pid = pid;
	localIDTable[atomId[i]].index = i;
      }
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
 *	$Author: brunner $	$Locker:  $		$State: Exp $
 *	$Revision: 1.1 $	$Date: 1996/10/24 18:50:33 $
 *
 ***************************************************************************
 * DESCRIPTION:
 *
 ***************************************************************************
 * REVISION HISTORY:
 *
 * $Log: AtomMap.C,v $
 * Revision 1.1  1996/10/24 18:50:33  brunner
 * Initial revision
 *
 *
 ***************************************************************************/

