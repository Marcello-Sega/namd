//-*-c++-*-
#ifndef COLLECTIONMASTER_H
#define COLLECTIONMASTER_H

#include "ckdefs.h"
#include "chare.h"
#include "c++interface.h"
#include "main.h"
#include "NamdTypes.h"
#include "ProcessorPrivate.h"

class CollectVectorMsg;

class CollectionMaster : public chare_object
{
public:

  static CollectionMaster *Object() { 
    return CpvAccess(CollectionMaster_instance); 
  }
  CollectionMaster(InitMsg *msg);
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
      seq(s), remaining(CNumPes()) { ; }

    // true -> send it and delete it!
    int append(AtomIDList &a, ResizeArray<Vector> &d)
    {
      int size = a.size();
      for( int i = 0; i < size; ++i )
      {
	data.add(VectorData(a[i],d[i]));
      }
      return ( ! --remaining );
    }

    int ready(void) { return ( ! remaining ); }

    int seq;

    operator<(const CollectVectorInstance &o) { return (seq < o.seq); }
    operator==(const CollectVectorInstance &o) { return (seq == o.seq); }
    void * operator new(size_t size) { return ::operator new(size); }
    void * operator new(size_t, void * ptr) { return ptr; }
    void operator delete(void* ptr) { ::operator delete(ptr); }

    class VectorData
    {
    public:
      VectorData(void) : aid(-1) { ; }
      VectorData(AtomID a, Vector d) : aid(a), data(d) { ; }
      AtomID aid;
      Vector data;
      operator<(const VectorData &o) { return (aid < o.aid); }
      operator==(const VectorData &o) { return (aid == o.aid); }
      void * operator new(size_t size) { return ::operator new(size); }
      void * operator new(size_t, void * ptr) { return ptr; }
      void operator delete(void* ptr) { ::operator delete(ptr); }
    };

    ResizeArray<VectorData> data;

  private:
    int remaining;

  };

  class CollectVectorSequence
  {
  public:

    CollectVectorInstance* submitData(
	int seq, AtomIDList &i, ResizeArray<Vector> &d)
    {
      CollectVectorInstance *c = data.find(CollectVectorInstance(seq));
      if ( ! c )
      {
	data.add(CollectVectorInstance(seq));
	c = data.find(CollectVectorInstance(seq));
      }
      if ( c->append(i,d) && queue.size() && queue[0] == seq )
      {
	c = new CollectVectorInstance(*c);
	data.del(CollectVectorInstance(seq));
	queue.del(0,1);
	return c;
      }
      else
      {
        return 0;
      }
    }

    CollectVectorInstance* enqueue(int seq)
    {
      queue.add(seq);
      if ( queue[0] == seq )
      {
        CollectVectorInstance *c = data.find(CollectVectorInstance(seq));
        if ( c && c->ready() )
        {
	  c = new CollectVectorInstance(*c);
	  data.del(CollectVectorInstance(seq));
	  queue.del(0,1);
	  return c;
        }
      }
      return 0;
    }

    ResizeArray<CollectVectorInstance> data;
    ResizeArray<int> queue;

    void * operator new(size_t size) { return ::operator new(size); }
    void operator delete(void* ptr) { ::operator delete(ptr); }
  };
private:

  CollectVectorSequence positions;
  CollectVectorSequence velocities;
  CollectVectorSequence forces;

};


class CollectVectorMsg : public comm_object
{
public:

  int seq;
  AtomIDList aid;
  ResizeArray<Vector> data;

  void * pack(int *length);
  void unpack(void *in);

  void * operator new(size_t s, int i, int j) {return comm_object::operator new(s,i,j);}
  void * operator new(size_t s, int i) {return comm_object::operator new(s,i);}
  void * operator new(size_t s) { return comm_object::operator new(s); }
  void * operator new(size_t, void *ptr) { return ptr; }

};

#endif

/***************************************************************************
 * RCS INFORMATION:
 *
 *	$RCSfile $
 *	$Author $	$Locker:  $		$State: Exp $
 *	$Revision: 1.1008 $	$Date: 1997/11/07 20:17:34 $
 *
 ***************************************************************************
 * REVISION HISTORY:
 *
 * $Log: CollectionMaster.h,v $
 * Revision 1.1008  1997/11/07 20:17:34  milind
 * Made NAMD to run on shared memory machines.
 *
 * Revision 1.1007  1997/07/09 21:26:39  milind
 * Ported NAMD2 to SP3. The SP specific code is within #ifdef SP2
 * and #endif's.
 *
 * Revision 1.1006  1997/04/06 22:44:55  ari
 * Add priorities to messages.  Mods to help proxies without computes.
 * Added quick enhancement to end of list insertion of ResizeArray(s)
 *
 * Revision 1.1005  1997/03/19 11:53:57  ari
 * Add Broadcast mechanism.
 * Fixed RCS Log entries on files that did not have Log entries.
 * Added some register variables to Molecule and ComputeNonbondedExcl.C
 *
 ***************************************************************************/
