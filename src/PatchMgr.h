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

class PatchMgrInitMsg;
class HomePatch;

class MovePatchesMsg : public comm_object {
public:
    NodeID  fromNodeID;
    PatchID pid;
    AtomIDList aid;
    PositionList p;
    VelocityList v;

    MovePatchesMsg(PatchID p, AtomIDList a, PositionList pl, VelocityList vl) : 
      pid(p), aid(a), p(pl), v(vl) {
      fromNodeID = CMyPe();
    }
    ~MovePatchesMsg() {};

    pack();
    unpack();
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

class HomePatch;

class HomePatchElem {
public:
  PatchID   pid;
  HomePatch *p;

  operator<(HomePatchElem e) { return (pid < e.pid); }
  operator==(HomePatchElem e) { return (pid == e.pid); }

  HomePatchElem(PatchID id=-1, HomePatch *patch=NULL) : pid(id), p(patch) {};
  ~HomePatchElem() { };
  HomePatchElem& operator=(const HomePatchElem &e) { 
    pid = e.pid; p = e.p;  // Do not delete p!  This op only used to shuffle
			   // we delete the p here only when the HomePatch is 
			   // moved off!
    return(*this);
  };
};

typedef SortedArray<HomePatchElem> HomePatchList;

struct MovePatch 
{
    MovePatch(NodeID n=-1, PatchID p=-1) : nodeID(n), pid(p) {};
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

class PatchMgr : public groupmember
{
private:
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
  
  int workDistribGroup;

public:
  PatchMgr(PatchMgrInitMsg *);
  ~PatchMgr();
  
  void recvMovePatches(MovePatchesMsg *msg);

  void createHomePatch(PatchID pid, AtomIDList aid, PositionList p, 
     VelocityList v); 
  void movePatch(PatchID, NodeID);
  void sigWorkDistrib();
  void sendMovePatches();
  void ackMovePatches(AckMovePatchesMsg *msg);
  HomePatch *homePatch(PatchID pid) {
     return homePatches.find(HomePatchElem(pid))->p;
  }
};

  // this is a list of immediate neighbor nodes (that is neighbor nodes
  // that has an immediate neighbor patch. We will wait for an
  // atom redistribution message from each node in this list per cycle.
  // And also we send an atom redistribution message to each of these nodes
  // even if no atoms moves there. Basically, these nodes must be in sync
  // as far as atom redist is concerned.
  // IntList   *neighborNodes;
  // void build_neighborlist();

  // check if atoms moved to another patch
  // handle the ones moved internally,
  // if an atom moved to a remote patch, send a msg to the node
  // this function is
  // void redistribute_atom();

  // receive messages from other nodes about atoms moving into patches 
  // void recv_atom_redist(AtomRedistMsg *msg);
  // each patch invokes this method to initiate atom redistribution
  // at the end of the cyle
  // void redist_my_atoms(int patchid);

#endif /* PATCHMGR_H */
/***************************************************************************
 * RCS INFORMATION:
 *
 *	$RCSfile: PatchMgr.h,v $
 *	$Author: ari $	$Locker:  $		$State: Exp $
 *	$Revision: 1.4 $	$Date: 1996/11/01 21:20:45 $
 *
 ***************************************************************************
 * REVISION HISTORY:
 *
 * $Log: PatchMgr.h,v $
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
 *
 ***************************************************************************/
