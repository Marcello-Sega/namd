/**
***  Copyright (c) 1995, 1996, 1997, 1998, 1999, 2000 by
***  The Board of Trustees of the University of Illinois.
***  All rights reserved.
**/

#include "CollectionMgr.decl.h"
#include "CollectionMgr.h"
#include "CollectionMaster.decl.h"
#include "CollectionMaster.h"
#include "Node.h"
#include "Molecule.h"
#include "SimParameters.h"

//#define DEBUGM
#include "Debug.h"

CollectionMgr::CollectionMgr(SlaveInitMsg *msg) : master(msg->master)
{
  delete msg;
  if (CpvAccess(CollectionMgr_instance) == 0) {
    CpvAccess(CollectionMgr_instance) = this;
  } else {
    DebugM(1, "CollectionMgr::CollectionMgr() - another instance of CollectionMgr exists!\n");
  }
}


CollectionMgr::~CollectionMgr(void)
{
}


void CollectionMgr::submitPositions(int seq, FullAtomList &a,
				Lattice l, int prec)
{
  Molecule *mol = Node::Object()->molecule;
  int numAtoms = a.size();
  AtomIDList aid(numAtoms);
  PositionList d(numAtoms);
  for ( int i=0; i<numAtoms; ++i ) {
    aid[i] = a[i].id;
    d[i] = l.reverse_transform(a[i].position,a[i].transform);
  }
  CollectVectorInstance *c;
  if ( ( c = positions.submitData(seq,aid,d,prec) ) )
  {
    CollectVectorMsg * msg = new CollectVectorMsg;
    msg->seq = c->seq;
    msg->aid = c->aid;
    msg->data = c->data;
    msg->fdata = c->fdata;
    CProxy_CollectionMaster cm(master);
    cm.receivePositions(msg);
    c->free();
  }
}


void CollectionMgr::submitVelocities(int seq, int zero, FullAtomList &a)
{
  int numAtoms = a.size();
  AtomIDList aid(numAtoms);
  PositionList d(numAtoms);
  for ( int i=0; i<numAtoms; ++i ) {
    aid[i] = a[i].id;
    if ( zero ) d[i] = 0.;
    else d[i] = a[i].velocity;
  }
  CollectVectorInstance *c;
  if ( ( c = velocities.submitData(seq,aid,d) ) )
  {
    CollectVectorMsg * msg = new CollectVectorMsg;
    msg->seq = c->seq;
    msg->aid = c->aid;
    msg->data = c->data;
    CProxy_CollectionMaster cm(master);
    cm.receiveVelocities(msg);
    c->free();
  }
}


void CollectionMgr::sendDataStream(const char *data) {
  DataStreamMsg *msg = new DataStreamMsg;
  msg->data.resize(strlen(data)+1);
  strcpy(msg->data.begin(),data);
  CProxy_CollectionMaster cm(master);
  cm.receiveDataStream(msg);
}


#include "CollectionMgr.def.h"

