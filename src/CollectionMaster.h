/**
***  Copyright (c) 1995, 1996, 1997, 1998, 1999, 2000 by
***  The Board of Trustees of the University of Illinois.
***  All rights reserved.
**/

#ifndef COLLECTIONMASTER_H
#define COLLECTIONMASTER_H

#include "charm++.h"
#include "main.h"
#include "NamdTypes.h"
#include "Lattice.h"
#include "ProcessorPrivate.h"
#include "PatchMap.h"
#include "PatchMap.inl"
#include "CollectionMaster.decl.h"
#include <stdio.h>

class CollectVectorMsg;
class DataStreamMsg;
class EnqueueDataMsg;

class CollectionMaster : public Chare
{
public:

  static CollectionMaster *Object() { 
    return CkpvAccess(CollectionMaster_instance); 
  }
  CollectionMaster();
  ~CollectionMaster(void);

  void receivePositions(CollectVectorMsg *msg);
  void receiveVelocities(CollectVectorMsg *msg);

  void receiveDataStream(DataStreamMsg *msg);

  void enqueuePositions(int seq, Lattice &lattice);
  void enqueueVelocities(int seq);

  void enqueuePositionsFromHandler(EnqueueDataMsg *msg);
  void enqueueVelocitiesFromHandler(int seq);

  class CollectVectorInstance;
  void disposePositions(CollectVectorInstance *c);
  void disposeVelocities(CollectVectorInstance *c);

  class CollectVectorInstance
  {
  public:

    CollectVectorInstance(void) : seq(-10) { ; }

    CollectVectorInstance(int s) { reset(s); }

    void free() { seq = -10; }
    int notfree() { return ( seq != -10 ); }

    void reset(int s) {
	if ( s == -10 ) NAMD_bug("seq == free in CollectionMaster");
        seq = s;
        remaining = (PatchMap::Object())->numNodesWithPatches(); 
	data.resize(0);
	fdata.resize(0);
    }

    // true -> send it and delete it!
    void append(AtomIDList &a, ResizeArray<Vector> &d, ResizeArray<FloatVector> &fd)
    {
      int size = a.size();
      if ( d.size() ) {
	for( int i = 0; i < size; ++i ) { data.item(a[i]) = d[i]; }
      }
      if ( fd.size() ) {
	for( int i = 0; i < size; ++i ) { fdata.item(a[i]) = fd[i]; }
      }
      --remaining;
    }

    int ready(void) { return ( ! remaining ); }

    int seq;
    Lattice lattice;

    ResizeArray<Vector> data;
    ResizeArray<FloatVector> fdata;

  private:
    int remaining;

  };

  class CollectVectorSequence
  {
  public:

    void submitData(
	int seq, AtomIDList &i, ResizeArray<Vector> &d, ResizeArray<FloatVector> &fd)
    {
      CollectVectorInstance **c = data.begin();
      CollectVectorInstance **c_e = data.end();
      for( ; c != c_e && (*c)->seq != seq; ++c );
      if ( c == c_e )
      {
        c = data.begin();
        for( ; c != c_e && (*c)->notfree(); ++c );
        if ( c == c_e ) {
          data.add(new CollectVectorInstance(seq));
          c = data.end() - 1;
        }
        (*c)->reset(seq);
      }
      (*c)->append(i,d,fd);
    }

    void enqueue(int seq, Lattice &lattice) {
      queue.add(seq);
      latqueue.add(lattice);
    }

    CollectVectorInstance* removeReady(void)
    {
      CollectVectorInstance *o = 0;
      if ( queue.size() )
      {
        int seq = queue[0];
        CollectVectorInstance **c = data.begin();
        CollectVectorInstance **c_e = data.end();
        for( ; c != c_e && (*c)->seq != seq; ++c );
        if ( c != c_e && (*c)->ready() )
        {
	  o = *c;
	  o->lattice = latqueue[0];
	  queue.del(0,1);
	  latqueue.del(0,1);
        }
      }
      return o;
    }

    ResizeArray<CollectVectorInstance*> data;
    ResizeArray<int> queue;
    ResizeArray<Lattice> latqueue;

  };
private:

  CollectVectorSequence positions;
  CollectVectorSequence velocities;
  FILE *dataStreamFile;

};


/* 
 * This class is used as a proxy to the CollectionMaster classi.
 * In the memory optimized version, the instance of CollectionMaster is 
 * moved to other processors. However, this instance is assumed to be on pe 0
 * and the Controller object contains the pointer to the instance. As a
 * result, a proxy class is needed to bridge the gap between the Controller
 * object on pe0 and the CollectionMaster object (which handles outputing
 * coordinates/velocities) on the other processor.
 */
class CollectionMasterHandler : public Chare
{
public:

  static CollectionMasterHandler *Object() { 
    return CkpvAccess(CollectionMasterHandler_instance); 
  }
  CollectionMasterHandler(MasterHandlerInitMsg *);
  ~CollectionMasterHandler(void);

  void enqueuePositions(EnqueueDataMsg *msg);
  void enqueueVelocities(int seq);
  void enqueuePositions(CkQdMsg *msg);
  void enqueueVelocities(CkQdMsg *msg);

  void setRealMaster(CkChareID m) { realMaster = m ;}
private:
  CkChareID realMaster;
  int enqueuePhase;
};

class CollectVectorMsg : public CMessage_CollectVectorMsg
{
public:

  int seq;
  AtomIDList aid;
  ResizeArray<Vector> data;
  ResizeArray<FloatVector> fdata;

  static void* pack(CollectVectorMsg* msg);
  static CollectVectorMsg* unpack(void *ptr);

};


class DataStreamMsg : public CMessage_DataStreamMsg {
public:

  ResizeArray<char> data;

  static void* pack(DataStreamMsg* msg);
  static DataStreamMsg* unpack(void *ptr);

};

class EnqueueDataMsg: public CMessage_EnqueueDataMsg {
public:
  int timestep;
  Lattice l;
  static void* pack(EnqueueDataMsg* msg);
  static EnqueueDataMsg* unpack(void *ptr);
};

class MasterHandlerInitMsg: public CMessage_MasterHandlerInitMsg {
public:
    CkChareID master;
};

#endif

