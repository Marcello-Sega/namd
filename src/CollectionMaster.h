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
    int append(AtomIDList &a, ResizeArray<Vector> &d)
    {
      int size = a.size();
      for( int i = 0; i < size; ++i )
      {
	data.item(a[i]) = d[i];
      }
      return ( ! --remaining );
    }

    int ready(void) { return ( ! remaining ); }

    int seq;

    int operator<(const CollectVectorInstance &o) { return (seq < o.seq); }
    int operator==(const CollectVectorInstance &o) { return (seq == o.seq); }
    void * operator new(size_t size) { return ::operator new(size); }
    void * operator new(size_t, void * ptr) { return ptr; }
    void operator delete(void* ptr) { ::operator delete(ptr); }
    // crashes egcs!  void operator delete(void*, void*) { ; }

    ResizeArray<Vector> data;

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


class CollectVectorMsg : public CMessage_CollectVectorMsg
{
public:

  int seq;
  AtomIDList aid;
  ResizeArray<Vector> data;

  static void* pack(CollectVectorMsg* msg);
  static CollectVectorMsg* unpack(void *ptr);

};

#endif

/***************************************************************************
 * RCS INFORMATION:
 *
 *	$RCSfile $
 *	$Author $	$Locker:  $		$State: Exp $
 *	$Revision: 1.1014 $	$Date: 1999/05/11 23:56:15 $
 *
 ***************************************************************************
 * REVISION HISTORY:
 *
 * $Log: CollectionMaster.h,v $
 * Revision 1.1014  1999/05/11 23:56:15  brunner
 * Changes for new charm version
 *
 * Revision 1.1013  1999/03/17 17:59:22  jim
 * Eliminated compiler warnings and errors.
 *
 * Revision 1.1012  1998/11/30 04:12:34  krishnan
 * Fixed the numNodes > nPatches bug.
 *
 * Revision 1.1011  1998/09/14 16:11:34  jim
 * Changes to reduce node 0 memory use.  Fixed bug in ResizeArray::item().
 *
 * Revision 1.1010  1998/03/03 23:05:02  brunner
 * Changed include files for new simplified Charm++ include file structure.
 *
 * Revision 1.1009  1998/01/15 05:40:53  jim
 * Added int return type to comparison operators.
 *
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
