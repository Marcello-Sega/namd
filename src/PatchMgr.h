/**
***  Copyright (c) 1995, 1996, 1997, 1998, 1999, 2000 by
***  The Board of Trustees of the University of Illinois.
***  All rights reserved.
**/

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

