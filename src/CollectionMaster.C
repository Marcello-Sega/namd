/**
***  Copyright (c) 1995, 1996, 1997, 1998, 1999, 2000 by
***  The Board of Trustees of the University of Illinois.
***  All rights reserved.
**/

#include "charm++.h"
#include "CollectionMaster.decl.h"
#include "CollectionMaster.h"
#include "InfoStream.h"
#include "Node.h"
#include "Output.h"
#include "ProcessorPrivate.h"
#include "packmsg.h"
#include <stdio.h>

// #define DEBUGM
#include "Debug.h"

CollectionMaster::CollectionMaster()
{
  if (CpvAccess(CollectionMaster_instance) == 0) {
    CpvAccess(CollectionMaster_instance) = this;
  } else {
    DebugM(1, "CollectionMaster::CollectionMaster() - another instance of CollectionMaster exists!\n");
  }
  dataStreamFile = 0;
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
    c->free();
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
    c->free();
}


void CollectionMaster::receiveDataStream(DataStreamMsg *msg) {
    CkPrintf(msg->data.begin());
    if ( ! dataStreamFile ) {
      dataStreamFile = fopen("aux_data.txt","w");
      if ( ! dataStreamFile ) NAMD_die("Can't open data stream file!");
    }
    fprintf(dataStreamFile,"%s",msg->data.begin());
    fflush(dataStreamFile);
    delete msg;
}


PACK_MSG(CollectVectorMsg,
  PACK(seq);
  PACK_RESIZE(aid);
  PACK_RESIZE(data);
  PACK_RESIZE(fdata);
)

PACK_MSG(DataStreamMsg,
  PACK_RESIZE(data);
)


#include "CollectionMaster.def.h"

