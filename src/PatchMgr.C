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

static char ident[] = "@(#)$Header: /home/cvs/namd/cvsroot/namd2/src/PatchMgr.C,v 1.2 1996/08/29 00:50:42 ari Exp $";


#include "ckdefs.h"
#include "chare.h"
#include "c++interface.h"

#include "PatchMgr.top.h"
#include "PatchMgr.h"

#include "NamdTypes.h"
#include "Compute.h"
#include "HomePatch.h"

#include "main.top.h"
#include "main.h"

#include "WorkDistrib.top.h"
#include "WorkDistrib.h"




PatchMgr::PatchMgr(PatchMgrInitMsg *msg)
{
    workDistribGroup = msg->workDistribGroup;
    delete msg;
}

PatchMgr::~PatchMgr()
{
}

void PatchMgr::createHomePatch(PatchID pid, AtomIDList aid, 
	PositionList p, VelocityList v) 
{
    homePatches.load(HomePatchElem(pid,new HomePatch(pid, aid, p, v)));
}

void PatchMgr::movePatch(PatchID pid, NodeID nodeID) 
{
    move.load(MovePatch(pid,nodeID));
}

void PatchMgr::sendMovePatches() 
{
    ackMovePending = move.size();
    if (ackMovePending == 0) {   // tell local WorkDistrib we are
	sigWorkDistrib();        // done with patch moves
	return;
    }
    MovePatchListIter m(move);
    for ( m = m.begin(); m != m.end(); m++) {
      HomePatch *p = homePatch((*m).pid);
      MovePatchesMsg *msg = new (MsgIndex(MovePatchesMsg))
	MovePatchesMsg((*m).pid, p->atomIDList, p->p, p->v);

      // Sending to PatchMgr::recvMovePatches on remote node
      CSendMsgBranch(PatchMgr, recvMovePatches, msg, thisgroup, (*m).nodeID);

      // Deleting the HomePatchElem will call a destructor for clean up
      // but the msg elements are safe since they use a container template
      // that uses ref counting.
      homePatches.del(HomePatchElem((*m).pid)); 
    }
    move.resize(0);
}

void PatchMgr::recvMovePatches(MovePatchesMsg *msg) {
    // Tell sending PatchMgr we received MovePatchMsg
    AckMovePatchesMsg *ackmsg = 
      new (MsgIndex(AckMovePatchesMsg)) AckMovePatchesMsg;
    CSendMsgBranch(PatchMgr, ackMovePatches, ackmsg, thisgroup, msg->fromNodeID);

    // Make a new HomePatch
    createHomePatch(msg->pid, msg->aid, msg->p, msg->v);
    delete msg;
}
    

void PatchMgr::ackMovePatches(AckMovePatchesMsg *msg)
{
    delete msg;
    if (! --ackMovePending) 
	sigWorkDistrib();
}

void PatchMgr::sigWorkDistrib()
{
    // Send msg to WorkDistrib that all patchMoves are completed
    MovePatchDoneMsg *msg = new (MsgIndex(MovePatchDoneMsg)) MovePatchDoneMsg;
    CSendMsgBranch(WorkDistrib, movePatchDone, msg, workDistribGroup, CMyPe());
}
   

#include "PatchMgr.bot.h"

/***************************************************************************
 * RCS INFORMATION:
 *
 *	$RCSfile: PatchMgr.C,v $
 *	$Author: ari $	$Locker:  $		$State: Exp $
 *	$Revision: 1.2 $	$Date: 1996/08/29 00:50:42 $
 *
 * REVISION HISTORY:
 *
 * $Log: PatchMgr.C,v $
 * Revision 1.2  1996/08/29 00:50:42  ari
 * *** empty log message ***
 *
 * Revision 1.1  1996/08/19 22:07:49  ari
 * Initial revision
 *
 ***************************************************************************/
