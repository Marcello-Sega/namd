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

class InitMsg;
class HomePatch;

typedef ResizeArray<HomePatch *> HomePatchList;

// This class is used to send atom redistribution information.
// It is derived from the Charm++ communication object
class AtomRedistMsg : public comm_object {
public:
    int x;
   
    AtomRedistMsg();
    ~AtomRedistMsg();
    pack();
    unpack();
};

typedef ResizeArray<int> LocalIndex;
typedef ResizeArray<Patch *> homePatchList;

// PatchMgr creates and manages homepatches. There exist one instance of 
// PatchMgr on each node (derived from Charm++ groupmember).
// That is, when a new operator causes creation of one instance on each node. 
// In addition to creation of homepatches, it handles the atom redistribution
// at the end of each cycle (i.e., atoms can move from patch to patch at the
// cycle boundaries).
class PatchMgr : public groupmember
{

private:
  int numAllPatches;
  int numHomePatches;

  // global patch number to local patch table conversion table
  LocalIndex indexGlb;

  // an array of patch pointers residing on this node
  HomePatchList homePatch;

public:
  PatchMgr(InitMsg *);
  ~PatchMgr();
  
  void createPatch(PatchID, PositionList&, VelocityList&);
  void movePatch(PatchIDList&, NodeID);
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
 *	$Revision: 1.1 $	$Date: 1996/08/19 22:07:49 $
 *
 ***************************************************************************
 * REVISION HISTORY:
 *
 * $Log: PatchMgr.h,v $
 * Revision 1.1  1996/08/19 22:07:49  ari
 * Initial revision
 *
 *
 ***************************************************************************/
