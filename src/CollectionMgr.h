//-*-c++-*-
#ifndef COLLECTIONMGR_H
#define COLLECTIONMGR_H

#include "ckdefs.h"
#include "chare.h"
#include "c++interface.h"
#include "main.h"
#include "NamdTypes.h"
#include "BOCgroup.h"
#include "PatchMap.h"


class CollectionMgr : public BOCclass
{
public:

  static CollectionMgr *Object() { return _instance; }
  CollectionMgr(SlaveInitMsg *msg);
  ~CollectionMgr(void);

  void submitPositions(int seq, AtomIDList &i, PositionList &d);
  void submitVelocities(int seq, AtomIDList &i, VelocityList &d);
  void submitForces(int seq, AtomIDList &i, ForceList &d);

private:

  static CollectionMgr *_instance;

  ChareIDType master;

  class CollectVectorInstance
  {
  public:

    CollectVectorInstance(void) : seq(-1) { ; }

    CollectVectorInstance(int s) :
      seq(s), remaining(PatchMap::Object()->numHomePatches()) { ; }

    // true -> send it and delete it!
    int append(AtomIDList &a, ResizeArray<Vector> &d)
    {
      int size = a.size();
      for( int i = 0; i < size; ++i )
      {
	aid.add(a[i]);
	data.add(d[i]);
      }
      return ( ! --remaining );
    }

    int seq;
    AtomIDList aid;
    ResizeArray<Vector> data;

    operator<(const CollectVectorInstance &o) { return (seq < o.seq); }
    operator==(const CollectVectorInstance &o) { return (seq == o.seq); }
    void * operator new(size_t size) { return ::operator new(size); }
    void operator delete(void* ptr) { ::operator delete(ptr); }

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
      if ( c->append(i,d) )
      {
	c = new CollectVectorInstance(*c);
	data.del(CollectVectorInstance(seq));
        return c;
      }
      else
      {
        return 0;
      }
    }

    ResizeArray<CollectVectorInstance> data;

    void * operator new(size_t size) { return ::operator new(size); }
    void operator delete(void* ptr) { ::operator delete(ptr); }
  };

  CollectVectorSequence positions;
  CollectVectorSequence velocities;
  CollectVectorSequence forces;

};


#endif
