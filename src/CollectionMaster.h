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
#include "ProcessorPrivate.h"
#include "PatchMap.h"
#include "PatchMap.inl"
#include "CollectionMaster.decl.h"
#include <stdio.h>

class CollectVectorMsg;
class DataStreamMsg;

class CollectionMaster : public Chare
{
public:

  static CollectionMaster *Object() { 
    return CpvAccess(CollectionMaster_instance); 
  }
  CollectionMaster();
  ~CollectionMaster(void);

  void receivePositions(CollectVectorMsg *msg);
  void receiveVelocities(CollectVectorMsg *msg);

  void receiveDataStream(DataStreamMsg *msg);

  void enqueuePositions(int seq);
  void enqueueVelocities(int seq);

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
        remaining = CkNumPes();
	int npatches=(PatchMap::Object())->numPatches(); 
	if (remaining > npatches) remaining=npatches; 
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

    void enqueue(int seq) { queue.add(seq); }

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
	  queue.del(0,1);
        }
      }
      return o;
    }

    ResizeArray<CollectVectorInstance*> data;
    ResizeArray<int> queue;

  };
private:

  CollectVectorSequence positions;
  CollectVectorSequence velocities;
  FILE *dataStreamFile;

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

#endif

