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

#ifndef PATCHMGR_H
#define PATCHMGR_H

#include "charm++.h"

#include "NamdTypes.h"
#include "SortedArray.h"
#include "HomePatch.h"
#include "HomePatchList.h"
#include "BOCgroup.h"
#include "Migration.h"
#include "MigrateAtomsMsg.h"
#include "PatchMgr.decl.h"


class HomePatch;

class MovePatchesMsg : public CMessage_MovePatchesMsg {
public:
    NodeID  fromNodeID;
    PatchID pid;
    AtomIDList aid;
    TransformList t;
    PositionList p;
    VelocityList v;

    MovePatchesMsg(void) { ; }

    MovePatchesMsg(PatchID n, AtomIDList a, TransformList tl,
				PositionList pl, VelocityList vl) : 
      pid(n), aid(a), t(tl), p(pl), v(vl)
    {
      fromNodeID = CkMyPe();
    }

  // pack and unpack functions
  static void* pack(MovePatchesMsg *msg);
  static MovePatchesMsg* unpack(void *ptr);
};

class MoveAtomMsg : public CMessage_MoveAtomMsg {
public:
  int atomid;
  int moveto;
  Vector coord;
};

// PatchMgr creates and manages homepatches. There exist one instance of 
// PatchMgr on each node (derived from Charm++ Group).  // That is, when a new operator causes creation of one instance on each node. 
// In addition to creation of homepatches, it handles the atom redistribution
// at the end of each cycle (i.e., atoms can move from patch to patch at the
// cycle boundaries).
struct MovePatch 
{
    MovePatch(PatchID p=-1, NodeID n=-1) : nodeID(n), pid(p) {};
    ~MovePatch() {};

    NodeID nodeID;
    PatchID pid;

    int operator<(MovePatch m) {
      return ( nodeID < m.nodeID );
    }

    int operator==(MovePatch m) {
      return ( nodeID == m.nodeID );
    }
};

typedef SortedArray<MovePatch> MovePatchList;
typedef ResizeArrayIter<MovePatch> MovePatchListIter;

class PatchMgr : public BOCclass
{

public:
  PatchMgr();
  ~PatchMgr();

  static PatchMgr* Object() { return CpvAccess(PatchMgr_instance); }
  

  void createHomePatch(PatchID pid, AtomIDList aid, TransformList t,
     PositionList p, VelocityList v); 

  void movePatch(PatchID, NodeID);
  void sendMovePatches();
  void recvMovePatches(MovePatchesMsg *msg);

  // void ackMovePatches(AckMovePatchesMsg *msg);

  HomePatch *homePatch(PatchID pid) {
     return homePatches.find(HomePatchElem(pid))->patch;
  } 

  void sendMigrationMsg(PatchID, MigrationInfo);
  void sendMigrationMsgs(PatchID, MigrationInfo*, int);
  void recvMigrateAtoms(MigrateAtomsMsg *);
  void recvMigrateAtomsCombined(MigrateAtomsCombinedMsg *);
  static void setGroup(BOCgroup g);

  void moveAtom(MoveAtomMsg *msg);
 
private:

  friend class PatchMap;
  PatchMap *patchMap;

  int numAllPatches;
  int numHomePatches;

  // an array of patch pointers residing on this node
  HomePatchList homePatches;

  // an array of patches to move off this node
  MovePatchList move;
  int ackMovePending;

  // data for combining migration messages
  MigrateAtomsCombinedMsg ** combineMigrationMsgs;
  int migrationCountdown;
};



#endif /* PATCHMGR_H */
/***************************************************************************
 * RCS INFORMATION:
 *
 *	$RCSfile: PatchMgr.h,v $
 *	$Author: jim $	$Locker:  $		$State: Exp $
 *	$Revision: 1.1013 $	$Date: 1999/08/11 16:53:09 $
 *
 ***************************************************************************
 * REVISION HISTORY:
 *
 * $Log: PatchMgr.h,v $
 * Revision 1.1013  1999/08/11 16:53:09  jim
 * Added move command to TCL scripting.
 *
 * Revision 1.1012  1999/05/11 23:56:44  brunner
 * Changes for new charm version
 *
 * Revision 1.1011  1998/08/11 16:30:30  jim
 * Modified output from periodic boundary simulations to return atoms to
 * internally consistent coordinates.  We store the transformations which
 * were performed and undo them at the end.  It might be better to do this
 * by always keeping the original coordinates and only doing the transform
 * for the nonbonded terms but this works for now.
 *
 * Revision 1.1010  1998/03/03 23:05:23  brunner
 * Changed include files for new simplified Charm++ include file structure.
 *
 * Revision 1.1009  1998/01/15 04:58:49  jim
 * Corrected "friend foo" to "friend class foo".
 *
 * Revision 1.1008  1997/12/26 23:10:56  milind
 * Made namd2 to compile, link and run under linux. Merged Templates and src
 * directoriies, and removed separate definition and declaration files for
 * templates.
 *
 * Revision 1.1007  1997/11/07 20:17:45  milind
 * Made NAMD to run on shared memory machines.
 *
 * Revision 1.1006  1997/04/11 06:03:29  jim
 * Message combining implemented for atom migration.
 *
 * Revision 1.1005  1997/04/10 22:29:17  jim
 * First steps towards combining atom migration messages.
 *
 * Revision 1.1004  1997/03/06 22:06:09  ari
 * Removed Compute.ci
 * Comments added - more code cleaning
 *
 * Revision 1.1003  1997/02/26 16:53:15  ari
 * Cleaning and debuging for memory leaks.
 * Adding comments.
 * Removed some dead code due to use of Quiescense detection.
 *
 * Revision 1.1002  1997/02/17 23:47:04  ari
 * Added files for cleaning up atom migration code
 *
 * Revision 1.1001  1997/02/13 04:43:13  jim
 * Fixed initial hanging (bug in PatchMap, but it still shouldn't have
 * happened) and saved migration messages in the buffer from being
 * deleted, but migration still dies (even on one node).
 *
 * Revision 1.1000  1997/02/06 15:59:07  ari
 * Resetting CVS to merge branches back into the main trunk.
 * We will stick to main trunk development as suggested by CVS manual.
 * We will set up tags to track fixed points of development/release
 * as suggested by CVS manual - all praise the CVS manual.
 *
 * Revision 1.779  1997/02/06 15:53:24  ari
 * Updating Revision Line, getting rid of branches
 *
 * Revision 1.778.2.1  1997/02/05 22:18:21  ari
 * Added migration code - Currently the framework is
 * there with compiling code.  This version does
 * crash shortly after migration is complete.
 * Migration appears to complete, but Patches do
 * not appear to be left in a correct state.
 *
 * Revision 1.778  1997/01/28 00:31:14  ari
 * internal release uplevel to 1.778
 *
 * Revision 1.777.2.2  1997/01/27 22:45:34  ari
 * Basic Atom Migration Code added.
 * Added correct magic first line to .h files for xemacs to go to C++ mode.
 * Compiles and runs without migration turned on.
 *
 * Revision 1.777.2.1  1997/01/21 23:04:51  ari
 * Basic framework for atom migration placed into code.  - Non
 * functional since it is not called.  Works currently without
 * atom migration.
 *
 * Revision 1.777  1997/01/17 19:36:50  ari
 * Internal CVS leveling release.  Start development code work
 * at 1.777.1.1.
 *
 * Revision 1.8  1996/12/16 23:46:01  jim
 * added placement new and explicit destructor calls to message
 *
 * Revision 1.7  1996/12/13 08:56:04  jim
 * move pack and unpack into C file, eliminated need for constructor
 * before unpack or destructor after pack
 *
 * Revision 1.6  1996/12/12 21:03:15  jim
 * added pack and unpack for patch movement message
 *
 * Revision 1.5  1996/11/22 00:18:51  ari
 * *** empty log message ***
 *
 * Revision 1.4  1996/11/01 21:20:45  ari
 * *** empty log message ***
 *
 * Revision 1.3  1996/09/03 22:54:25  ari
 * *** empty log message ***
 *
 * Revision 1.2  1996/08/29 00:50:42  ari
 * *** empty log message ***
 *
 * Revision 1.1  1996/08/19 22:07:49  ari
 * Initial revision
 *
 ***************************************************************************/
