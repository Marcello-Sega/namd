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
  inline static AtomMap *Object() {return _instance;}
  ~AtomMap(void);

  void allocateMap(int nAtomIDs);

  int registerIDs(PatchID pid, AtomIDList al);

  LocalID localID(AtomID id);

  void clearMap(void);


protected:
  AtomMap(void);

private:
  static AtomMap *_instance;

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
 *	$Author: ari $	$Locker:  $		$State: Exp $
 *	$Revision: 1.777 $	$Date: 1997/01/17 19:35:18 $
 *
 ***************************************************************************
 * REVISION HISTORY:
 *
 * $Log: AtomMap.h,v $
 * Revision 1.777  1997/01/17 19:35:18  ari
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

