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

static char ident[] = "@(#)$Header: /home/cvs/namd/cvsroot/namd2/src/PatchMgr.C,v 1.1005 1997/03/06 22:06:08 ari Exp $";


#include "ckdefs.h"
#include "chare.h"
#include "c++interface.h"

#include "PatchMgr.top.h"
#include "PatchMgr.h"

#include "NamdTypes.h"
//#include "Compute.h"
#include "HomePatch.h"
#include "PatchMap.h"

#include "main.top.h"
#include "main.h"

#include "WorkDistrib.top.h"
#include "WorkDistrib.h"

// #include "Node.top.h"
// #include "Node.h"

// #define DEBUGM
#define MIN_DEBUG_LEVEL 3
#include "Debug.h"


// Singleton pattern - static initialization
PatchMgr *PatchMgr::_instance = 0;

// BOC constructor
PatchMgr::PatchMgr(InitMsg *msg)
{
    delete msg;

    // Singleton pattern
    if (_instance == NULL) {
	_instance = this;
    } else {
	iout << iFILE << iERROR << iPE 
	  << "PatchMgr instanced twice on same node!" << endi;
	CharmExit();
    }

    // Get PatchMap singleton started
    patchMap = PatchMap::Instance();
    patchMap->registerPatchMgr(this);
}

PatchMgr::~PatchMgr()
{
    HomePatchListIter hi(homePatches);
    for ( hi = hi.begin(); hi != hi.end(); hi++) {
      HomePatchElem* elem = homePatches.find(HomePatchElem(hi->pid));
      delete elem->patch;
    }
}


void PatchMgr::createHomePatch(PatchID pid, AtomIDList aid, 
	PositionList p, VelocityList v) 
{
    HomePatch *patch = new HomePatch(pid, aid, p, v);
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

      MovePatchesMsg *msg = new (MsgIndex(MovePatchesMsg))
	MovePatchesMsg(m->pid, p->atomIDList, p->p, p->v);

      // Sending to PatchMgr::recvMovePatches on remote node
      CSendMsgBranch(PatchMgr, recvMovePatches, msg, thisgroup, m->nodeID);

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
    createHomePatch(msg->pid, msg->aid, msg->p, msg->v);
    delete msg;

    // Tell sending PatchMgr we received MovePatchMsg
//    AckMovePatchesMsg *ackmsg = 
//      new (MsgIndex(AckMovePatchesMsg)) AckMovePatchesMsg;
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
  MigrateAtomsMsg *msg = 
    new (MsgIndex(MigrateAtomsMsg)) MigrateAtomsMsg(src,m.destPatchID,m.mList);
  CSendMsgBranch(PatchMgr,recvMigrateAtoms,msg,thisgroup,m.destNodeID);
}

// Receive end of sendMigrateionMsg() above
void PatchMgr::recvMigrateAtoms (MigrateAtomsMsg *msg) {
  //  msg must be deleted by HomePatch::depositMigrationMsg();
  PatchMap::Object()->homePatch(msg->destPatchID)->depositMigration(msg);
}

void * MovePatchesMsg::pack (int *length)
  {
    DebugM(1,"MovePatchesMsg::pack() - aid.size() = " << aid.size() << endl);
    DebugM(1,"MovePatchesMsg::pack() - p.size() = " << p.size() << endl);
    DebugM(1,"MovePatchesMsg::pack() - v.size() = " << v.size() << endl);
    *length = sizeof(NodeID) + sizeof(PatchID) + sizeof(int) +
		aid.size() * sizeof(AtomID) +
		p.size() * sizeof(Position) +
		v.size() * sizeof(Velocity);
    char *buffer = (char*)new_packbuffer(this,*length);
    char *b = buffer;
    *((NodeID*)b) = fromNodeID; b += sizeof(NodeID);
    *((PatchID*)b) = pid; b += sizeof(PatchID);
    *((int*)b) = aid.size(); b += sizeof(int);
    for ( int i = 0; i < aid.size(); i++ )
    {
      *((AtomID*)b) = aid[i]; b += sizeof(AtomID);
      *((Position*)b) = p[i]; b += sizeof(Position);
      *((Velocity*)b) = v[i]; b += sizeof(Velocity);
    }
    this->~MovePatchesMsg();
    return buffer;
  }

void MovePatchesMsg::unpack (void *in)
  {
    new((void*)this) MovePatchesMsg;
    char *b = (char*)in;
    fromNodeID = *((NodeID*)b); b += sizeof(NodeID);
    pid = *((PatchID*)b); b += sizeof(PatchID);
    int size = *((int*)b); b += sizeof(int);
    DebugM(1,"MovePatchesMsg::unpack() - size = " << size << endl);
    aid.resize(size);
    p.resize(size);
    v.resize(size);
    for ( int i = 0; i < size; i++ )
    {
      aid[i] = *((AtomID*)b); b += sizeof(AtomID);
      p[i] = *((Position*)b); b += sizeof(Position);
      v[i] = *((Velocity*)b); b += sizeof(Velocity);
    }
  }

#include "PatchMgr.bot.h"

/***************************************************************************
 * RCS INFORMATION:
 *
 *	$RCSfile: PatchMgr.C,v $
 *	$Author: ari $	$Locker:  $		$State: Exp $
 *	$Revision: 1.1005 $	$Date: 1997/03/06 22:06:08 $
 *
 * REVISION HISTORY:
 *
 * $Log: PatchMgr.C,v $
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
