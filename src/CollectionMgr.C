#include "CollectionMgr.top.h"
#include "CollectionMgr.h"
#include "CollectionMaster.top.h"
#include "CollectionMaster.h"


CollectionMgr::CollectionMgr(SlaveInitMsg *msg) : master(msg->master)
{
  delete msg;
}


CollectionMgr::~CollectionMgr(void)
{
}


void CollectionMgr::submitPositions(int seq, AtomIDList &i, PositionList &d)
{
  CollectVectorInstance *c;
  if ( c = positions.submitData(seq,i,d) )
  {
    CollectVectorMsg * msg = new (MsgIndex(CollectVectorMsg)) CollectVectorMsg;
    msg->seq = c->seq;
    msg->aid = c->aid;
    msg->data = c->data;
    CSendMsg(CollectionMaster,receivePositions,msg,&master);
    delete c;
  }
}


void CollectionMgr::submitVelocities(int seq, AtomIDList &i, VelocityList &d)
{
  CollectVectorInstance *c;
  if ( c = velocities.submitData(seq,i,d) )
  {
    CollectVectorMsg * msg = new (MsgIndex(CollectVectorMsg)) CollectVectorMsg;
    msg->seq = c->seq;
    msg->aid = c->aid;
    msg->data = c->data;
    CSendMsg(CollectionMaster,receiveVelocities,msg,&master);
    delete c;
  }
}


void CollectionMgr::submitForces(int seq, AtomIDList &i, ForceList &d)
{
  CollectVectorInstance *c;
  if ( c = forces.submitData(seq,i,d) )
  {
    CollectVectorMsg * msg = new (MsgIndex(CollectVectorMsg)) CollectVectorMsg;
    msg->seq = c->seq;
    msg->aid = c->aid;
    msg->data = c->data;
    CSendMsg(CollectionMaster,receiveForces,msg,&master);
    delete c;
  }
}


#include "CollectionMgr.bot.h"

