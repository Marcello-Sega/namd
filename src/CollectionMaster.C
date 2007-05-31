/**
***  Copyright (c) 1995, 1996, 1997, 1998, 1999, 2000 by
***  The Board of Trustees of the University of Illinois.
***  All rights reserved.
**/

#include "largefiles.h"

#include "InfoStream.h"
#include "CollectionMaster.h"
#include "Node.h"
#include "Output.h"
#include "ProcessorPrivate.h"
#include "SimParameters.h"
#include "packmsg.h"
#include "CollectionMaster.decl.h"

// #define DEBUGM
#include "Debug.h"

//CollectionMasterHandler should always be on processor 0
CollectionMasterHandler::CollectionMasterHandler(MasterHandlerInitMsg *msg): realMaster(msg->master)
{
  delete msg;
  if (CpvAccess(CollectionMasterHandler_instance) == 0) {
    CpvAccess(CollectionMasterHandler_instance) = this;
  } else {
    DebugM(1, "CollectionMasterHandler::CollectionMasterHandler() - another instance of CollectionMasterHandler exists!\n");
  }
  enqueuePhase = 0;
}


CollectionMasterHandler::~CollectionMasterHandler(void)
{
}

void CollectionMasterHandler::enqueuePositions(EnqueueDataMsg *msg){
    if(enqueuePhase==0){
	CProxy_CollectionMaster cm(realMaster);
	EnqueueDataMsg *newmsg = new EnqueueDataMsg;
	newmsg->timestep = msg->timestep;
	newmsg->l = msg->l;
	cm.enqueuePositionsFromHandler(newmsg);
	delete msg;
    }else if(enqueuePhase==1){
	enqueuePhase = 0;
	#if CHARM_VERSION>050402
	CkStartQD(CkIndex_CollectionMasterHandler::enqueuePositions((CkQdMsg*)0), &thishandle);
	#else
	CkStartQD(CProxy_CollectionMasterHandler::ckIdx_enqueuePositions((CkQdMsg*)0), &thishandle);
	#endif
    }else{
	NAMD_die("Enqueue phase at enqueuePositions in the CollectionMasterHandler has wrong value!\n");
    }
}

void CollectionMasterHandler::enqueueVelocities(int seq){
    if(enqueuePhase==0){
	CProxy_CollectionMaster cm(realMaster);
	cm.enqueueVelocitiesFromHandler(seq);
    }else if(enqueuePhase==1){
	enqueuePhase = 0;
	#if CHARM_VERSION>050402
	CkStartQD(CkIndex_CollectionMasterHandler::enqueueVelocities((CkQdMsg*)0), &thishandle);
	#else
	CkStartQD(CProxy_CollectionMasterHandler::ckIdx_enqueueVelocities((CkQdMsg*)0), &thishandle);
	#endif
    }else{
	NAMD_die("Enqueue phase at enqueueVelocities in the CollectionMasterHandler has wrong value!\n");
    }
}

void CollectionMasterHandler::enqueuePositions(CkQdMsg *qmsg){
    delete qmsg;
    Object()->enqueuePositions((EnqueueDataMsg *)NULL);
}

void CollectionMasterHandler::enqueueVelocities(CkQdMsg *qmsg){
    delete qmsg;
    Object()->enqueueVelocities(0);
}

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

void CollectionMaster::enqueuePositions(int seq, Lattice &lattice)
{
  positions.enqueue(seq,lattice);

  CollectVectorInstance *c;
  while ( ( c = positions.removeReady() ) ) { disposePositions(c); }
}

void CollectionMaster::enqueuePositionsFromHandler(EnqueueDataMsg *msg){
    enqueuePositions(msg->timestep, msg->l);
    delete msg;
}

void CollectionMaster::disposePositions(CollectVectorInstance *c)
{
    DebugM(3,"Collected positions at " << c->seq << std::endl);
    int seq = c->seq;
    int size = c->data.size();
    if ( ! size ) size = c->fdata.size();
    Vector *data = c->data.begin();
    FloatVector *fdata = c->fdata.begin();
    Node::Object()->output->coordinate(seq,size,data,fdata,c->lattice);
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
  Lattice dummy;
  velocities.enqueue(seq,dummy);

  CollectVectorInstance *c;
  while ( ( c = velocities.removeReady() ) ) { disposeVelocities(c); }
}

void CollectionMaster::enqueueVelocitiesFromHandler(int seq){
    enqueueVelocities(seq);
}

void CollectionMaster::disposeVelocities(CollectVectorInstance *c)
{
    DebugM(3,"Collected velocities at " << c->seq << std::endl);
    int seq = c->seq;
    int size = c->data.size();
    Vector *data = c->data.begin();
    Node::Object()->output->velocity(seq,size,data);
    c->free();
}


void CollectionMaster::receiveDataStream(DataStreamMsg *msg) {
    if ( ! dataStreamFile ) {
      char *fname = Node::Object()->simParameters->auxFilename;
      // iout has large file linking issues on AIX
      // iout << iINFO << "OPENING AUXILIARY DATA STREAM FILE "
      // 				<< fname << "\n" << endi;
      CkPrintf("Info: OPENING AUXILIARY DATA STREAM FILE %s\n", fname);
      NAMD_backup_file(fname);
      dataStreamFile = fopen(fname,"w");
      if ( ! dataStreamFile )
		NAMD_die("Can't open auxiliary data stream file!");
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

PACK_MSG(EnqueueDataMsg,
  PACK(timestep);
  PACK(l);
)

#include "CollectionMaster.def.h"

