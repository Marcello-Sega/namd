/***************************************************************************/
/*                                                                         */
/*              (C) Copyright 1996 The Board of Trustees of the            */
/*                          University of Illinois                         */
/*                           All Rights Reserved                           */
/*                                                                         */
/***************************************************************************/

/***************************************************************************/
/* DESCRIPTION:                                                            */
/*								           */
/***************************************************************************/

#include "charm++.h"

#include "PatchMgr.decl.h"
#include "PatchMgr.h"

#include "NamdTypes.h"
//#include "Compute.h"
#include "HomePatch.h"
#include "PatchMap.h"

#include "main.decl.h"
#include "main.h"

#include "WorkDistrib.decl.h"
#include "WorkDistrib.h"

// #include "Node.decl.h"
// #include "Node.h"

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
  // We note that m.mList may be NULL indicating no atoms to migrate
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

void * MovePatchesMsg::pack (MovePatchesMsg *m)
  {
    DebugM(1,"MovePatchesMsg::pack() - aid.size() = " << m->aid.size() << endl);
    DebugM(1,"MovePatchesMsg::pack() - p.size() = " << m->p.size() << endl);
    DebugM(1,"MovePatchesMsg::pack() - v.size() = " << m->v.size() << endl);
    int length = sizeof(NodeID) + sizeof(PatchID) + sizeof(int) +
		m->aid.size() * sizeof(AtomID) +
		m->t.size() * sizeof(Transform) +
		m->p.size() * sizeof(Position) +
		m->v.size() * sizeof(Velocity);
    char *buffer = (char*)CkAllocBuffer(m,length);
    char *b = buffer;
    memcpy(b, &(m->fromNodeID), sizeof(NodeID)); b += sizeof(NodeID);
    memcpy(b, &(m->pid), sizeof(PatchID)); b += sizeof(PatchID);
    int size=m->aid.size(); 
    memcpy(b, &size, sizeof(int)); b += sizeof(int);
    memcpy(b, m->aid.begin(),size*sizeof(AtomID)); b += size*sizeof(AtomID);
    memcpy(b, m->t.begin(),size*sizeof(Transform)); b += size*sizeof(Transform);
    memcpy(b, m->p.begin(),size*sizeof(Position)); b += size*sizeof(Position);
    memcpy(b, m->v.begin(),size*sizeof(Velocity)); b += size*sizeof(Velocity);
    delete m;
    return buffer;
  }

MovePatchesMsg* MovePatchesMsg::unpack (void *ptr)
  {
    void *_ptr = CkAllocBuffer(ptr, sizeof(MovePatchesMsg));
    MovePatchesMsg* m = new (_ptr) MovePatchesMsg;
    char *b = (char*)ptr;
    memcpy(&(m->fromNodeID), b, sizeof(NodeID)); b += sizeof(NodeID);
    memcpy(&(m->pid),b, sizeof(PatchID)); b += sizeof(PatchID);
    int size; memcpy(&size, b, sizeof(int)); b += sizeof(int);
    DebugM(1,"MovePatchesMsg::unpack() - size = " << size << endl);
    m->aid.resize(size);
    memcpy(m->aid.begin(),b,size*sizeof(AtomID)); b += size*sizeof(AtomID);
    m->t.resize(size);
    memcpy(m->t.begin(),b,size*sizeof(Transform)); b += size*sizeof(Transform);
    m->p.resize(size);
    memcpy(m->p.begin(),b,size*sizeof(Position)); b += size*sizeof(Position);
    m->v.resize(size);
    memcpy(m->v.begin(),b,size*sizeof(Velocity)); b += size*sizeof(Velocity);
    CkFreeMsg(ptr);
    return m;
  }

#include "PatchMgr.def.h"

/***************************************************************************
 * RCS INFORMATION:
 *
 *	$RCSfile: PatchMgr.C,v $
 *	$Author: brunner $	$Locker:  $		$State: Exp $
 *	$Revision: 1.1017 $	$Date: 1999/05/11 23:56:43 $
 *
 * REVISION HISTORY:
 *
 * $Log: PatchMgr.C,v $
 * Revision 1.1017  1999/05/11 23:56:43  brunner
 * Changes for new charm version
 *
 * Revision 1.1016  1998/10/24 19:57:52  jim
 * Eliminated warnings generated by g++ -Wall.
 *
 * Revision 1.1015  1998/09/13 21:06:12  jim
 * Cleaned up output, defaults, etc.
 *
 * Revision 1.1014  1998/08/11 16:30:30  jim
 * Modified output from periodic boundary simulations to return atoms to
 * internally consistent coordinates.  We store the transformations which
 * were performed and undo them at the end.  It might be better to do this
 * by always keeping the original coordinates and only doing the transform
 * for the nonbonded terms but this works for now.
 *
 * Revision 1.1013  1998/03/03 23:05:22  brunner
 * Changed include files for new simplified Charm++ include file structure.
 *
 * Revision 1.1012  1998/02/10 23:30:30  milind
 * Fixed to reflect the current changes to Charm++ translator.
 *
 * Revision 1.1011  1997/12/19 23:42:37  jim
 * Replaced assignments with memcpys and reordered memcpys for efficiency.
 *
 * Revision 1.1010  1997/11/07 20:17:45  milind
 * Made NAMD to run on shared memory machines.
 *
 * Revision 1.1009  1997/04/11 06:03:26  jim
 * Message combining implemented for atom migration.
 *
 * Revision 1.1008  1997/04/10 22:29:16  jim
 * First steps towards combining atom migration messages.
 *
 * Revision 1.1007  1997/04/10 09:14:03  ari
 * Final debugging for compute migration / proxy creation for load balancing.
 * Lots of debug code added, mostly turned off now.
 * Fixed bug in PositionBox when Patch had no dependencies.
 * Eliminated use of cout and misuse of iout in numerous places.
 *                                            Ari & Jim
 *
 * Revision 1.1006  1997/04/04 23:34:25  milind
 * Got NAMD2 to run on Origin2000.
 * Included definitions of class static variables in C files.
 * Fixed alignment bugs by using memcpy instead of assignment in
 * pack and unpack.
 *
 * Revision 1.1005  1997/03/06 22:06:08  ari
 * Removed Compute.ci
 * Comments added - more code cleaning
 *
 * Revision 1.1004  1997/02/26 16:53:14  ari
 * Cleaning and debuging for memory leaks.
 * Adding comments.
 * Removed some dead code due to use of Quiescense detection.
 *
 * Revision 1.1003  1997/02/17 23:47:03  ari
 * Added files for cleaning up atom migration code
 *
 * Revision 1.1002  1997/02/13 04:43:12  jim
 * Fixed initial hanging (bug in PatchMap, but it still shouldn't have
 * happened) and saved migration messages in the buffer from being
 * deleted, but migration still dies (even on one node).
 *
 * Revision 1.1001  1997/02/07 07:50:22  jim
 * Removed debug messages.
 *
 * Revision 1.1000  1997/02/06 15:59:05  ari
 * Resetting CVS to merge branches back into the main trunk.
 * We will stick to main trunk development as suggested by CVS manual.
 * We will set up tags to track fixed points of development/release
 * as suggested by CVS manual - all praise the CVS manual.
 *
 * Revision 1.779  1997/02/06 15:53:23  ari
 * Updating Revision Line, getting rid of branches
 *
 * Revision 1.778.2.1  1997/02/05 22:18:20  ari
 * Added migration code - Currently the framework is
 * there with compiling code.  This version does
 * crash shortly after migration is complete.
 * Migration appears to complete, but Patches do
 * not appear to be left in a correct state.
 *
 * Revision 1.778  1997/01/28 00:31:12  ari
 * internal release uplevel to 1.778
 *
 * Revision 1.777.2.1  1997/01/21 23:04:49  ari
 * Basic framework for atom migration placed into code.  - Non
 * functional since it is not called.  Works currently without
 * atom migration.
 *
 * Revision 1.777  1997/01/17 19:36:48  ari
 * Internal CVS leveling release.  Start development code work
 * at 1.777.1.1.
 *
 * Revision 1.9  1996/12/19 00:37:31  jim
 * increase MIN_DEBUG_LEVEL
 *
 * Revision 1.8  1996/12/16 23:46:01  jim
 * added placement new and explicit destructor calls to message
 *
 * Revision 1.7  1996/12/13 08:56:04  jim
 * move pack and unpack into C file, eliminated need for constructor
 * before unpack or destructor after pack
 *
 * Revision 1.6  1996/11/22 00:18:51  ari
 * *** empty log message ***
 *
 * Revision 1.5  1996/11/04 17:13:13  ari
 * *** empty log message ***
 *
 * Revision 1.4  1996/11/01 21:20:45  ari
 * *** empty log message ***
 *
 * Revision 1.3  1996/09/03 22:51:14  ari
 * *** empty log message ***
 *
 * Revision 1.2  1996/08/29 00:50:42  ari
 * *** empty log message ***
 *
 * Revision 1.1  1996/08/19 22:07:49  ari
 * Initial revision
 *
 ***************************************************************************/
