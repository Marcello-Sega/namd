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
  inline static AtomMap *Object() {return _instance;}
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
  static AtomMap *_instance;

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

/***************************************************************************
 * RCS INFORMATION:
 *
 *	$RCSfile: AtomMap.h,v $
 *	$Author: ari $	$Locker:  $		$State: Exp $
 *	$Revision: 1.1004 $	$Date: 1997/04/10 09:13:47 $
 *
 ***************************************************************************
 * REVISION HISTORY:
 *
 * $Log: AtomMap.h,v $
 * Revision 1.1004  1997/04/10 09:13:47  ari
 * Final debugging for compute migration / proxy creation for load balancing.
 * Lots of debug code added, mostly turned off now.
 * Fixed bug in PositionBox when Patch had no dependencies.
 * Eliminated use of cout and misuse of iout in numerous places.
 *                                            Ari & Jim
 *
 * Revision 1.1003  1997/03/04 22:37:03  ari
 * Clean up of code.  Debug statements removal, dead code removal.
 * Minor fixes, output fixes.
 * Commented some code from the top->down.  Mainly reworked Namd, Node, main.
 *
 * Revision 1.1002  1997/02/13 16:17:10  ari
 * Intermediate debuging commit - working to fix deep bug in migration?
 *
 * Revision 1.1001  1997/02/07 05:42:29  ari
 * Some bug fixing - atom migration on one node works
 * Atom migration on multiple nodes gets SIGSEGV
 *
 * Revision 1.1000  1997/02/06 15:57:28  ari
 * Resetting CVS to merge branches back into the main trunk.
 * We will stick to main trunk development as suggested by CVS manual.
 * We will set up tags to track fixed points of development/release
 * as suggested by CVS manual - all praise the CVS manual.
 *
 * Revision 1.778  1997/01/28 00:29:47  ari
 * internal release uplevel to 1.778
 *
 * Revision 1.777.2.1  1997/01/27 22:44:49  ari
 * Basic Atom Migration Code added.
 * Added correct magic first line to .h files for xemacs to go to C++ mode.
 * Compiles and runs without migration turned on.
 *
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

