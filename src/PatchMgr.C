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

static char ident[] = "@(#)$Header: /home/cvs/namd/cvsroot/namd2/src/PatchMgr.C,v 1.779 1997/02/06 15:53:23 ari Exp $";


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
#define MIN_DEBUG_LEVEL 4
#include "Debug.h"


PatchMgr::PatchMgr(InitMsg *msg)
{
    delete msg;

    _instance = this;
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
    DebugM(2,"sendMovePatches() - moving " << move.size() << " patches.\n");
    ackMovePending = move.size();
    if (ackMovePending == 0) {   // tell local WorkDistrib we are
	WorkDistrib::messageMovePatchDone();  // done with patch moves
	return;
    }
    MovePatchListIter m(move);
    for ( m = m.begin(); m != m.end(); m++) {
      DebugM(1,"sendMovePatches() - moving patch " << (*m).pid << ".\n");
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

     DebugM(1,"received patch " << msg->pid << " from node " << msg->fromNodeID << ".\n");
     DebugM(1,"recvMovePatches() - creating patch " << msg->pid << ".\n");
     DebugM(1,"aid.size() = " << msg->aid.size() << "\n");
     DebugM(1,"p.size() = " << msg->p.size() << "\n");
     DebugM(1,"v.size() = " << msg->v.size() << "\n");
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

void PatchMgr::recvMigrateAtoms (MigrateAtomsMsg *msg) {
  DebugM(4, "Received Migration Msg from node " << msg->fromNodeID << "\n");
  DebugM(4, "         Migration Msg from patch " << msg->srcPatchID << "\n");
  DebugM(4, "         Migration Msg to patch " << msg->destPatchID << "\n");
  PatchMap::Object()->homePatch(msg->destPatchID)->depositMigration(msg->srcPatchID,
    msg->migrationList);
  delete msg;
}

void PatchMgr::sendMigrationMsg(PatchID src, MigrationInfo m) {
  DebugM(3, "Received Migration List from " << src << " size = " << m->size() << "\n" );
  
  MigrateAtomsMsg *msg = new (MsgIndex(MigrateAtomsMsg)) 
     MigrateAtomsMsg(src,m.destPatchID,m.mList);
  if (m.mList) {
    DebugM(4,"Sending "<<m.destPatchID<<" from "<<src<<"size="<<(*m.mList).size()<<"\n" );
  } else {
    DebugM(4,"Sending "<<m.destPatchID<<" from "<<src<<" nothing\n");
  }
  CSendMsgBranch(PatchMgr,recvMigrateAtoms,msg,thisgroup,m.destNodeID);
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



void * MigrateAtomsMsg::pack (int *length) {
    if (migrationList != NULL) {
      DebugM(4,"MigrateAtomsMsg::pack() - migrationList->size() = " << migrationList->size() << "\n" );
      *length = sizeof(NodeID) + sizeof(PatchID) + sizeof(PatchID)
	      + sizeof(int) + migrationList->size() * sizeof(MigrationElem);
    } else {
      DebugM(4,"MigrateAtomsMsg::pack() NULL message\n" );
      *length = sizeof(NodeID) + sizeof(PatchID) + sizeof(PatchID)
	      +	sizeof(int);
    }
    char *buffer = (char*)new_packbuffer(this,*length);
    char *b = buffer;
    *((NodeID*)b) = fromNodeID; b += sizeof(NodeID);
    *((PatchID*)b) = srcPatchID; b += sizeof(PatchID);
    *((PatchID*)b) = destPatchID; b += sizeof(PatchID);

    if (migrationList != NULL) {
      *((int*)b) = migrationList->size(); b += sizeof(int);
      for ( int i = 0; i < migrationList->size(); i++ )
      {
	*((MigrationElem*)b) = (*migrationList)[i]; b += sizeof(MigrationElem);
      }
    }
    else {
      *((int*)b) = 0; b += sizeof(int);
    }

    this->~MigrateAtomsMsg();
    return buffer;
}

void MigrateAtomsMsg::unpack (void *in) {
  new((void*)this) MigrateAtomsMsg;
  char *b = (char*)in;
  fromNodeID = *((NodeID*)b); b += sizeof(NodeID);
  srcPatchID = *((PatchID*)b); b += sizeof(PatchID);
  destPatchID = *((PatchID*)b); b += sizeof(PatchID);
  int size = *((int*)b); b += sizeof(int);
  DebugM(4,"MigrateAtomsMsg::unpack() - from node = " << fromNodeID << endl);
  DebugM(4,"MigrateAtomsMsg::unpack() - from patch = " << srcPatchID << endl);
  DebugM(4,"MigrateAtomsMsg::unpack() - to patch = " << destPatchID << endl);
  DebugM(4,"MigrateAtomsMsg::unpack() - size = " << size << endl);
  if (size != 0) {
    migrationList = new MigrationList();
    migrationList->resize(size);
    for ( int i = 0; i < size; i++ )
    {
      (*migrationList)[i] = *((MigrationElem*)b); b += sizeof(MigrationElem);
    }
  }
  else {
    migrationList = NULL;
  }
}


#include "PatchMgr.bot.h"

/***************************************************************************
 * RCS INFORMATION:
 *
 *	$RCSfile: PatchMgr.C,v $
 *	$Author: ari $	$Locker:  $		$State: Exp $
 *	$Revision: 1.779 $	$Date: 1997/02/06 15:53:23 $
 *
 * REVISION HISTORY:
 *
 * $Log: PatchMgr.C,v $
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
