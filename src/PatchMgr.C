/**
***  Copyright (c) 1995, 1996, 1997, 1998, 1999, 2000 by
***  The Board of Trustees of the University of Illinois.
***  All rights reserved.
**/

#include "charm++.h"

#include "PatchMgr.decl.h"
#include "PatchMgr.h"

#include "NamdTypes.h"
//#include "Compute.h"
#include "HomePatch.h"
#include "PatchMap.h"
#include "AtomMap.h"

#include "main.decl.h"
#include "main.h"

#include "WorkDistrib.decl.h"
#include "WorkDistrib.h"

#include "packmsg.h"

// #define DEBUGM
#define MIN_DEBUG_LEVEL 3
#include "Debug.h"


// BOC constructor
PatchMgr::PatchMgr()
{
    // CkPrintf("[%d] PatchMgr Created\n", CkMyPe());

    // Singleton pattern
    if (CpvAccess(PatchMgr_instance) == NULL) {
	CpvAccess(PatchMgr_instance) = this;
    } else {
	iout << iFILE << iERROR << iPE 
	  << "PatchMgr instanced twice on same processor!" << endi;
	CkExit();
    }

    // Get PatchMap singleton started
    patchMap = PatchMap::Instance();
    patchMap->registerPatchMgr(this);

    // Message combining initialization
    migrationCountdown = 0;
    combineMigrationMsgs = new MigrateAtomsCombinedMsg*[CkNumPes()];
}

PatchMgr::~PatchMgr()
{
    HomePatchListIter hi(homePatches);
    for ( hi = hi.begin(); hi != hi.end(); hi++) {
      HomePatchElem* elem = homePatches.find(HomePatchElem(hi->pid));
      delete elem->patch;
    }
    delete [] combineMigrationMsgs;
}


void PatchMgr::createHomePatch(PatchID pid, AtomIDList aid, 
	TransformList t, PositionList p, VelocityList v) 
{
    HomePatch *patch = new HomePatch(pid, aid, t, p, v);
    homePatches.load(HomePatchElem(pid, patch));
    patchMap->registerPatch(pid, patch);
}


// Add a HomePatch to a list of patches to be moved 
// HomePatches are actually moved by invoking sendMovePatches() below
void PatchMgr::movePatch(PatchID pid, NodeID nodeID) 
{
    move.load(MovePatch(pid,nodeID));
}


// Uses list constructed by movePatch() and dispatches
// HomePatch(es) to new nodes
void PatchMgr::sendMovePatches() 
{
    if (! move.size())
	return;

    MovePatchListIter m(move);
    for ( m = m.begin(); m != m.end(); m++) {
      HomePatch *p = homePatch(m->pid);
      patchMap->unregisterPatch(m->pid, p);

      MovePatchesMsg *msg = new 
	MovePatchesMsg(m->pid, p->atomIDList, p->t, p->p, p->v);

      // Sending to PatchMgr::recvMovePatches on remote node
      CProxy_PatchMgr cp(thisgroup);
      cp.recvMovePatches(msg, m->nodeID);

      // Deleting the HomePatchElem will call a destructor for clean up
      // but the msg elements are safe since they use a container template
      // that uses ref counting.
      delete p;
      homePatches.del(HomePatchElem(m->pid)); 
    }
    move.resize(0);
}

void PatchMgr::recvMovePatches(MovePatchesMsg *msg) {
    // Make a new HomePatch
    createHomePatch(msg->pid, msg->aid, msg->t, msg->p, msg->v);
    delete msg;

    // Tell sending PatchMgr we received MovePatchMsg
//    AckMovePatchesMsg *ackmsg = 
//      new AckMovePatchesMsg;
//    CSendMsgBranch(PatchMgr,ackMovePatches, ackmsg, thisgroup, msg->fromNodeID);
}
    

//void PatchMgr::ackMovePatches(AckMovePatchesMsg *msg)
//{
//    delete msg;
//    if (! --ackMovePending) 
//	WorkDistrib::messageMovePatchDone();
//}



// Called by HomePatch to migrate atoms off to new patches
// Message combining could occur here
void PatchMgr::sendMigrationMsg(PatchID src, MigrationInfo m) {
  MigrateAtomsMsg *msg = new MigrateAtomsMsg(src,m.destPatchID,m.mList);
  CProxy_PatchMgr cp(thisgroup);
  cp.recvMigrateAtoms(msg, m.destNodeID);
}

// Called by HomePatch to migrate atoms off to new patches
// Message combining occurs here
void PatchMgr::sendMigrationMsgs(PatchID src, MigrationInfo *m, int numMsgs) {
/*
  for (int i=0; i < numMsgs; i++) {
    PatchMgr::Object()->sendMigrationMsg(src, m[i]);
  }
*/
  if ( ! migrationCountdown )  // (re)initialize
  {
    // DebugM(3,"migrationCountdown (re)initialize\n");
    numHomePatches = patchMap->numHomePatches();
    migrationCountdown = numHomePatches;
    int numPes = CkNumPes();
    for ( int i = 0; i < numPes; ++i ) combineMigrationMsgs[i] = 0;
  }
  for (int i=0; i < numMsgs; i++) {  // buffer messages
    int destNodeID = m[i].destNodeID;
    if ( 1 ) // destNodeID != CkMyPe() )
    {
      if ( ! combineMigrationMsgs[destNodeID] )
      {
        combineMigrationMsgs[destNodeID] = new MigrateAtomsCombinedMsg();
      }
      combineMigrationMsgs[destNodeID]->add(src,m[i].destPatchID,m[i].mList);
    }
    else
    {
	// for now buffer local messages too
    }
  }
  migrationCountdown -= 1;
  // DebugM(3,"migrationCountdown = " << migrationCountdown << "\n");
  if ( ! migrationCountdown )  // send out combined messages
  {
    int numPes = CkNumPes();
    for ( int destNodeID = 0; destNodeID < numPes; ++destNodeID )
      if ( combineMigrationMsgs[destNodeID] )
      {
	DebugM(3,"Sending MigrateAtomsCombinedMsg to node " << destNodeID << "\n");
        CProxy_PatchMgr cp(thisgroup);
        cp.recvMigrateAtomsCombined(combineMigrationMsgs[destNodeID],destNodeID);
      }
  }
}

// Receive end of sendMigrationMsg() above
void PatchMgr::recvMigrateAtoms (MigrateAtomsMsg *msg) {
  //  msg must be deleted by HomePatch::depositMigrationMsg();
  PatchMap::Object()->homePatch(msg->destPatchID)->depositMigration(msg);
}

void PatchMgr::recvMigrateAtomsCombined (MigrateAtomsCombinedMsg *msg)
{
  DebugM(3,"Received MigrateAtomsCombinedMsg with " << msg->srcPatchID.size() << " messages.\n");
  msg->distribute();
  delete msg;
}

void PatchMgr::moveAtom(MoveAtomMsg *msg) {
  LocalID lid = AtomMap::Object()->localID(msg->atomid);
  if ( lid.pid != notUsed ) {
    HomePatch *hp = patchMap->homePatch(lid.pid);
    if ( hp ) {
      if ( msg->moveto ) { hp->p[lid.index] = msg->coord; }
      else { hp->p[lid.index] += msg->coord; }
    }
  }
  delete msg;
}


PACK_MSG(MovePatchesMsg,
  PACK(fromNodeID);
  PACK(pid);
  PACK_RESIZE(aid);
  PACK_RESIZE(t);
  PACK_RESIZE(p);
  PACK_RESIZE(v);
)


#include "PatchMgr.def.h"

