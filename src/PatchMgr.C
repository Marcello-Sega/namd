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

static char ident[] = "@(#)$Header: /home/cvs/namd/cvsroot/namd2/src/PatchMgr.C,v 1.6 1996/11/22 00:18:51 ari Exp $";


#include "ckdefs.h"
#include "chare.h"
#include "c++interface.h"

#include "PatchMgr.top.h"
#include "PatchMgr.h"

#include "NamdTypes.h"
#include "Compute.h"
#include "HomePatch.h"
#include "PatchMap.h"

#include "main.top.h"
#include "main.h"

#include "WorkDistrib.top.h"
#include "WorkDistrib.h"

#include "Node.top.h"
#include "Node.h"

#define DEBUGM
#include "Debug.h"


PatchMgr::PatchMgr(InitMsg *msg)
{
    delete msg;

    patchMap = PatchMap::Instance();
    patchMap->registerPatchMgr(this);
}

PatchMgr::~PatchMgr()
{
    // This is should clean up and delete the only thing that dangles,
    // the home patch pointers.  Of course, we expect PatchMgr to dissappear!
    // This is a kludge - the lists should have better semantics.

    MovePatchListIter m(move);
    for ( m = m.begin(); m != m.end(); m++) {
      HomePatchElem* hp = homePatches.find(HomePatchElem((*m).pid));
      delete hp->p;
    }
}


void PatchMgr::createHomePatch(PatchID pid, AtomIDList aid, 
	PositionList p, VelocityList v) 
{
    DebugM(3, "PatchMgr::createHomePatch() number atoms = " << aid.size() << endl );
    HomePatch *patch = new HomePatch(pid, aid, p, v);
    homePatches.load(HomePatchElem(pid, patch));
    patchMap->registerPatch(pid, patch);
}

void PatchMgr::movePatch(PatchID pid, NodeID nodeID) 
{
    move.load(MovePatch(pid,nodeID));
}

void PatchMgr::sendMovePatches() 
{
    ackMovePending = move.size();
    if (ackMovePending == 0) {   // tell local WorkDistrib we are
	WorkDistrib::messageMovePatchDone();  // done with patch moves
	return;
    }
    MovePatchListIter m(move);
    for ( m = m.begin(); m != m.end(); m++) {
      HomePatch *p = homePatch((*m).pid);
      patchMap->unregisterPatch((*m).pid, p);

      MovePatchesMsg *msg = new (MsgIndex(MovePatchesMsg))
	MovePatchesMsg((*m).pid, p->atomIDList, p->p, p->v);

      // Sending to PatchMgr::recvMovePatches on remote node
      CSendMsgBranch(PatchMgr, recvMovePatches, msg, thisgroup, (*m).nodeID);

      // Deleting the HomePatchElem will call a destructor for clean up
      // but the msg elements are safe since they use a container template
      // that uses ref counting.
      delete p;
      homePatches.del(HomePatchElem((*m).pid)); 
    }
    move.resize(0);
}

void PatchMgr::recvMovePatches(MovePatchesMsg *msg) {
    // Tell sending PatchMgr we received MovePatchMsg
    AckMovePatchesMsg *ackmsg = 
      new (MsgIndex(AckMovePatchesMsg)) AckMovePatchesMsg;
    CSendMsgBranch(PatchMgr,ackMovePatches, ackmsg, thisgroup, msg->fromNodeID);

    // Make a new HomePatch
    createHomePatch(msg->pid, msg->aid, msg->p, msg->v);
    delete msg;
}
    

void PatchMgr::ackMovePatches(AckMovePatchesMsg *msg)
{
    delete msg;
    if (! --ackMovePending) 
	WorkDistrib::messageMovePatchDone();
}


#include "PatchMgr.bot.h"

/***************************************************************************
 * RCS INFORMATION:
 *
 *	$RCSfile: PatchMgr.C,v $
 *	$Author: ari $	$Locker:  $		$State: Exp $
 *	$Revision: 1.6 $	$Date: 1996/11/22 00:18:51 $
 *
 * REVISION HISTORY:
 *
 * $Log: PatchMgr.C,v $
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
