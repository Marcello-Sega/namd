//-*-c++-*-
#ifndef COLLECTIONMASTER_H
#define COLLECTIONMASTER_H

#include "charm++.h"
#include "main.h"
#include "NamdTypes.h"
#include "ProcessorPrivate.h"
#include "PatchMap.h"
#include "PatchMap.inl"
#include "CollectionMaster.decl.h"

class CollectVectorMsg;

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
  void receiveForces(CollectVectorMsg *msg);

  void enqueuePositions(int seq);
  void enqueueVelocities(int seq);
  void enqueueForces(int seq);

  class CollectVectorInstance;
  void disposePositions(CollectVectorInstance *c);
  void disposeVelocities(CollectVectorInstance *c);
  void disposeForces(CollectVectorInstance *c);

  class CollectVectorInstance
  {
  public:

    CollectVectorInstance(void) : seq(-1) { ; }

    CollectVectorInstance(int s) :
      seq(s), remaining(CkNumPes()) { 
		int npatches=(PatchMap::Object())->numPatches(); 
		if (CkNumPes() > npatches) remaining=npatches; 
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

    int operator<(const CollectVectorInstance &o) { return (seq < o.seq); }
    int operator==(const CollectVectorInstance &o) { return (seq == o.seq); }

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
      CollectVectorInstance *c = data.find(CollectVectorInstance(seq));
      if ( ! c )
      {
	data.add(CollectVectorInstance(seq));
	c = data.find(CollectVectorInstance(seq));
      }
      c->append(i,d,fd);
    }

    void enqueue(int seq) { queue.add(seq); }

    CollectVectorInstance* removeReady(void)
    {
      CollectVectorInstance *o = 0;
      if ( queue.size() )
      {
        int seq = queue[0];
        CollectVectorInstance *c = data.find(CollectVectorInstance(seq));
        if ( c && c->ready() )
        {
	  o = new CollectVectorInstance(*c);
	  data.del(*c);
	  queue.del(0,1);
        }
      }
      return o;
    }

    ResizeArray<CollectVectorInstance> data;
    ResizeArray<int> queue;

  };
private:

  CollectVectorSequence positions;
  CollectVectorSequence velocities;
  CollectVectorSequence forces;

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

#endif

