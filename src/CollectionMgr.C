#include "charm++.h"
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


void CollectionMgr::submitPositions(int seq, AtomIDList &i, PositionList &d,
				Lattice l, TransformList &t, int prec)
{
  int wrapWater = Node::Object()->simParameters->wrapWater;
  Molecule *mol = Node::Object()->molecule;
  PositionList d2(d.size());
  PositionList::iterator d_i,d_e,d2_i;
  AtomIDList::iterator a_i = i.begin();
  TransformList::iterator t_i = t.begin();
  d_i = d.begin();  d_e = d.end();  d2_i = d2.begin();
  for ( ; d_i != d_e ; ++d_i, ++d2_i, ++a_i, ++t_i ) {
    if ( wrapWater && mol->is_water(*a_i) ) {
      *d2_i = *d_i;
    } else {
      *d2_i = l.reverse_transform(*d_i,*t_i);
    }
  }
  CollectVectorInstance *c;
  if ( ( c = positions.submitData(seq,i,d2,prec) ) )
  {
    CollectVectorMsg * msg = new CollectVectorMsg;
    msg->seq = c->seq;
    msg->aid = c->aid;
    msg->data = c->data;
    msg->fdata = c->fdata;
    CProxy_CollectionMaster cm(master);
    cm.receivePositions(msg);
    delete c;
  }
}


void CollectionMgr::submitVelocities(int seq, AtomIDList &i, VelocityList &d)
{
  CollectVectorInstance *c;
  if ( ( c = velocities.submitData(seq,i,d) ) )
  {
    CollectVectorMsg * msg = new CollectVectorMsg;
    msg->seq = c->seq;
    msg->aid = c->aid;
    msg->data = c->data;
    CProxy_CollectionMaster cm(master);
    cm.receiveVelocities(msg);
    delete c;
  }
}


void CollectionMgr::submitForces(int seq, AtomIDList &i, ForceList &d)
{
  CollectVectorInstance *c;
  if ( ( c = forces.submitData(seq,i,d) ) )
  {
    CollectVectorMsg * msg = new CollectVectorMsg;
    msg->seq = c->seq;
    msg->aid = c->aid;
    msg->data = c->data;
    CProxy_CollectionMaster cm(master);
    cm.receiveForces(msg);
    delete c;
  }
}


#include "CollectionMgr.def.h"

