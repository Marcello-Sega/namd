#include "charm++.h"
#include "CollectionMaster.decl.h"
#include "CollectionMaster.h"
#include "InfoStream.h"
#include "Node.h"
#include "Output.h"
#include "ProcessorPrivate.h"
#include "packmsg.h"

// #define DEBUGM
#include "Debug.h"

CollectionMaster::CollectionMaster()
{
  if (CpvAccess(CollectionMaster_instance) == 0) {
    CpvAccess(CollectionMaster_instance) = this;
  } else {
    DebugM(1, "CollectionMaster::CollectionMaster() - another instance of CollectionMaster exists!\n");
  }
}


CollectionMaster::~CollectionMaster(void)
{
}


void CollectionMaster::receivePositions(CollectVectorMsg *msg)
{
  positions.submitData(msg->seq,msg->aid,msg->data,msg->fdata);
  delete msg;

  CollectVectorInstance *c;
  while ( ( c = positions.removeReady() ) ) { disposePositions(c); }
}

void CollectionMaster::enqueuePositions(int seq)
{
  positions.enqueue(seq);

  CollectVectorInstance *c;
  while ( ( c = positions.removeReady() ) ) { disposePositions(c); }
}

void CollectionMaster::disposePositions(CollectVectorInstance *c)
{
    DebugM(3,"Collected positions at " << c->seq << endl);
    int seq = c->seq;
    int size = c->data.size();
    if ( ! size ) size = c->fdata.size();
    Vector *data = c->data.begin();
    FloatVector *fdata = c->fdata.begin();
    Node::Object()->output->coordinate(seq,size,data,fdata);
    delete c;
}


void CollectionMaster::receiveVelocities(CollectVectorMsg *msg)
{
  velocities.submitData(msg->seq,msg->aid,msg->data,msg->fdata);
  delete msg;

  CollectVectorInstance *c;
  while ( ( c = velocities.removeReady() ) ) { disposeVelocities(c); }
}

void CollectionMaster::enqueueVelocities(int seq)
{
  velocities.enqueue(seq);

  CollectVectorInstance *c;
  while ( ( c = velocities.removeReady() ) ) { disposeVelocities(c); }
}

void CollectionMaster::disposeVelocities(CollectVectorInstance *c)
{
    DebugM(3,"Collected velocities at " << c->seq << endl);
    int seq = c->seq;
    int size = c->data.size();
    Vector *data = c->data.begin();
    Node::Object()->output->velocity(seq,size,data);
    delete c;
}


void CollectionMaster::receiveForces(CollectVectorMsg *msg)
{
  forces.submitData(msg->seq,msg->aid,msg->data,msg->fdata);
  delete msg;

  CollectVectorInstance *c;
  while ( ( c = forces.removeReady() ) ) { disposeForces(c); }
}

void CollectionMaster::enqueueForces(int seq)
{
  forces.enqueue(seq);

  CollectVectorInstance *c;
  while ( ( c = forces.removeReady() ) ) { disposeForces(c); }
}

void CollectionMaster::disposeForces(CollectVectorInstance *c)
{
    DebugM(3,"Collected forces at " << c->seq << endl);
    int seq = c->seq;
    int size = c->data.size();
    Vector *data = c->data.begin();
    Node::Object()->output->all_force(seq,size,data);
    delete c;
}


PACK_MSG(CollectVectorMsg,
  PACK(seq);
  PACK_RESIZE(aid);
  PACK_RESIZE(data);
  PACK_RESIZE(fdata);
)


#include "CollectionMaster.def.h"

