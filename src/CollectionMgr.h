/**
***  Copyright (c) 1995, 1996, 1997, 1998, 1999, 2000 by
***  The Board of Trustees of the University of Illinois.
***  All rights reserved.
**/

#ifndef COLLECTIONMGR_H
#define COLLECTIONMGR_H

#include "charm++.h"

#include "main.h"
#include "NamdTypes.h"
#include "BOCgroup.h"
#include "PatchMap.h"
#include "ProcessorPrivate.h"
#include "CollectionMgr.decl.h"


class SlaveInitMsg : public CMessage_SlaveInitMsg
{
public:
  CkChareID master;
};

class CollectionMgr : public BOCclass
{
public:

  static CollectionMgr *Object() { 
    return CpvAccess(CollectionMgr_instance); 
  }
  CollectionMgr(SlaveInitMsg *msg);
  ~CollectionMgr(void);

  void submitPositions(int seq, AtomIDList &i, PositionList &d,
				Lattice l, TransformList &t, int prec);
  void submitVelocities(int seq, AtomIDList &i, VelocityList &d);

  class CollectVectorInstance
  {
  public:

    CollectVectorInstance(void) : seq(-1) { ; }

    CollectVectorInstance(int s) : seq(s) { ; }

    CollectVectorInstance(int s, int p) : seq(s), precisions(p),
      remaining(PatchMap::Object()->numHomePatches()) { ; }

    // true -> send it and delete it!
    int append(AtomIDList &a, ResizeArray<Vector> &d)
    {
      int size = a.size();
      for( int i = 0; i < size; ++i )
      {
	aid.add(a[i]);
	if ( precisions & 2 ) data.add(d[i]);
	if ( precisions & 1 ) fdata.add(d[i]);
      }
      return ( ! --remaining );
    }

    int seq;
    AtomIDList aid;
    int precisions;
    ResizeArray<Vector> data;
    ResizeArray<FloatVector> fdata;

    int operator<(const CollectVectorInstance &o) { return (seq < o.seq); }
    int operator==(const CollectVectorInstance &o) { return (seq == o.seq); }

  private:
    int remaining;

  };

  class CollectVectorSequence
  {
  public:

    CollectVectorInstance* submitData(
	int seq, AtomIDList &i, ResizeArray<Vector> &d, int prec=2)
    {
      CollectVectorInstance *c = data.find(CollectVectorInstance(seq));
      if ( ! c )
      {
	data.add(CollectVectorInstance(seq,prec));
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

  };
private:

  CkChareID master;


  CollectVectorSequence positions;
  CollectVectorSequence velocities;

};

#endif

