/***************************************************************************/
/*                                                                         */
/*              (C) Copyright 1996 The Board of Trustees of the            */
/*                          University of Illinois                         */
/*                           All Rights Reserved                           */
/*									   */
/***************************************************************************/

#ifndef ATOMMAP_H
#define ATOMMAP_H

#include "NamdTypes.h"

enum { notUsed = -1 };

class AtomMap
{
public:
  static AtomMap *handle;

  AtomMap(void);

  ~AtomMap(void);

  void allocateMap(int nAtomIDs);

  int registerIDs(PatchID pid, AtomID *atomid, int nIds);

  LocalID localID(AtomID id);

  void clearMap(void);

private:
  LocalID *localIDTable;
  int tableSz;
  Boolean cleared;

};

//----------------------------------------------------------------------

inline LocalID AtomMap::localID(AtomID id)
{
  return localIDTable[id];
}

#endif /* ATOMMAP_H */

/***************************************************************************
 * RCS INFORMATION:
 *
 *	$RCSfile: AtomMap.h,v $
 *	$Author: brunner $	$Locker:  $		$State: Exp $
 *	$Revision: 1.1 $	$Date: 1996/10/24 18:50:33 $
 *
 ***************************************************************************
 * DESCRIPTION:
 *
 ***************************************************************************
 * REVISION HISTORY:
 *
 * $Log: AtomMap.h,v $
 * Revision 1.1  1996/10/24 18:50:33  brunner
 * Initial revision
 *
 *
 ***************************************************************************/

