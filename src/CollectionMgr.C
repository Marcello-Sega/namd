#include "CollectionMgr.top.h"
#include "CollectionMgr.h"
#include "CollectionMaster.top.h"
#include "CollectionMaster.h"
#include "Priorities.h"

#define DEBUGM
#include "Debug.h"

CollectionMgr *CollectionMgr::_instance = 0;

CollectionMgr::CollectionMgr(SlaveInitMsg *msg) : master(msg->master)
{
  delete msg;
  if (_instance == 0) {
    _instance = this;
  } else {
    DebugM(1, "CollectionMgr::CollectionMgr() - another instance of CollectionMgr exists!\n");
  }
}


CollectionMgr::~CollectionMgr(void)
{
}


void CollectionMgr::submitPositions(int seq, AtomIDList &i, PositionList &d)
{
  CollectVectorInstance *c;
  if ( c = positions.submitData(seq,i,d) )
  {
    CollectVectorMsg * msg = 
      new (MsgIndex(CollectVectorMsg),Priorities::numBits) CollectVectorMsg;
    msg->seq = c->seq;
    msg->aid = c->aid;
    msg->data = c->data;
    //*CPriorityPtr(msg) = Priorities::low;
    //CSetQueueing(msg, C_QUEUEING_IFIFO);
    CSendMsg(CollectionMaster,receivePositions,msg,&master);
    delete c;
  }
}


void CollectionMgr::submitVelocities(int seq, AtomIDList &i, VelocityList &d)
{
  CollectVectorInstance *c;
  if ( c = velocities.submitData(seq,i,d) )
  {
    CollectVectorMsg * msg 
      = new (MsgIndex(CollectVectorMsg),Priorities::numBits) CollectVectorMsg;
    msg->seq = c->seq;
    msg->aid = c->aid;
    msg->data = c->data;
    //*CPriorityPtr(msg) = Priorities::low;
    //CSetQueueing(msg, C_QUEUEING_IFIFO);
    CSendMsg(CollectionMaster,receiveVelocities,msg,&master);
    delete c;
  }
}


void CollectionMgr::submitForces(int seq, AtomIDList &i, ForceList &d)
{
  CollectVectorInstance *c;
  if ( c = forces.submitData(seq,i,d) )
  {
    CollectVectorMsg * msg 
      = new (MsgIndex(CollectVectorMsg), Priorities::numBits) CollectVectorMsg;
    msg->seq = c->seq;
    msg->aid = c->aid;
    msg->data = c->data;
    //*CPriorityPtr(msg) = Priorities::low;
    //CSetQueueing(msg, C_QUEUEING_IFIFO);
    CSendMsg(CollectionMaster,receiveForces,msg,&master);
    delete c;
  }
}


#include "CollectionMgr.bot.h"

/***************************************************************************
 * RCS INFORMATION:
 *
 *	$RCSfile $
 *	$Author $	$Locker:  $		$State: Exp $
 *	$Revision: 1.1003 $	$Date: 1997/04/08 07:08:07 $
 *
 ***************************************************************************
 * REVISION HISTORY:
 *
 * $Log: CollectionMgr.C,v $
 * Revision 1.1003  1997/04/08 07:08:07  ari
 * Modification for dynamic loadbalancing - moving computes
 * Still bug in new computes or usage of proxies/homepatches.
 * Works if ldbStrategy is none as before.
 *
 * Revision 1.1002  1997/04/06 22:44:56  ari
 * Add priorities to messages.  Mods to help proxies without computes.
 * Added quick enhancement to end of list insertion of ResizeArray(s)
 *
 * Revision 1.1001  1997/03/19 11:53:58  ari
 * Add Broadcast mechanism.
 * Fixed RCS Log entries on files that did not have Log entries.
 * Added some register variables to Molecule and ComputeNonbondedExcl.C
 *
 ***************************************************************************/
