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

#include "ckdefs.h"
#include "chare.h"
#include "c++interface.h"

#include "NamdTypes.h"
#include "Templates/SortedArray.h"
#include "HomePatch.h"
#include "HomePatchList.h"
#include "BOCgroup.h"
#include "Migration.h"
#include "MigrateAtomsMsg.h"


class InitMsg;
class HomePatch;

class MovePatchesMsg : public comm_object {
public:
    NodeID  fromNodeID;
    PatchID pid;
    AtomIDList aid;
    PositionList p;
    VelocityList v;

    MovePatchesMsg(void) { ; }

    MovePatchesMsg(PatchID n, AtomIDList a, PositionList pl, VelocityList vl) : 
      pid(n), aid(a), p(pl), v(vl)
    {
      fromNodeID = CMyPe();
    }

  void * operator new(size_t s, int i) {return comm_object::operator new(s,i);}
  void * operator new(size_t s) { return comm_object::operator new(s); }
  void * operator new(size_t, void *ptr) { return ptr; }

  // pack and unpack functions
  void * pack (int *length);
  void unpack (void *in);
};

class AckMovePatchesMsg : public comm_object {
public:
   int dummy;

   AckMovePatchesMsg() {};
   ~AckMovePatchesMsg() {};
};


// PatchMgr creates and manages homepatches. There exist one instance of 
// PatchMgr on each node (derived from Charm++ groupmember).
// That is, when a new operator causes creation of one instance on each node. 
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

typedef ResizeArray<int> PatchIndex;

class PatchMgr : public BOCclass
{

public:
  PatchMgr(InitMsg *);
  ~PatchMgr();

  static PatchMgr* Object() { return _instance; }
  
  void recvMovePatches(MovePatchesMsg *msg);

  void createHomePatch(PatchID pid, AtomIDList aid, PositionList p, 
     VelocityList v); 
  void movePatch(PatchID, NodeID);
  void sendMovePatches();
  void ackMovePatches(AckMovePatchesMsg *msg);
  HomePatch *homePatch(PatchID pid) {
     return homePatches.find(HomePatchElem(pid))->p;
  } 

  void sendMigrationMsg(PatchID, MigrationInfo);
  void recvMigrateAtoms(MigrateAtomsMsg *);
  static void setGroup(BOCgroup g);
 
private:
  static PatchMgr* _instance;

  friend PatchMap;
  PatchMap *patchMap;

  int numAllPatches;
  int numHomePatches;

  // global patch number to local patch table conversion table
  PatchIndex patchIndex;

  // an array of patch pointers residing on this node
  HomePatchList homePatches;

  // an array of patches to move off this node
  MovePatchList move;
  int ackMovePending;
};



#endif /* PATCHMGR_H */
/***************************************************************************
 * RCS INFORMATION:
 *
 *	$RCSfile: PatchMgr.h,v $
 *	$Author: ari $	$Locker:  $		$State: Exp $
 *	$Revision: 1.1002 $	$Date: 1997/02/17 23:47:04 $
 *
 ***************************************************************************
 * REVISION HISTORY:
 *
 * $Log: PatchMgr.h,v $
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
