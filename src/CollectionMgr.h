//-*-c++-*-
#ifndef COLLECTIONMGR_H
#define COLLECTIONMGR_H

#include "charm++.h"

#include "main.h"
#include "NamdTypes.h"
#include "BOCgroup.h"
#include "PatchMap.h"
#include "ProcessorPrivate.h"
#include "CollectionMgr.decl.h"


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
  void submitForces(int seq, AtomIDList &i, ForceList &d);

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
  CollectVectorSequence forces;

};

#endif

/***************************************************************************
 * RCS INFORMATION:
 *
 *	$RCSfile $
 *	$Author $	$Locker:  $		$State: Exp $
 *	$Revision: 1.1012 $	$Date: 1999/09/12 19:33:16 $
 *
 ***************************************************************************
 * REVISION HISTORY:
 *
 * $Log: CollectionMgr.h,v $
 * Revision 1.1012  1999/09/12 19:33:16  jim
 * Collections now use floats when possible.
 *
 * Revision 1.1011  1999/07/06 20:32:40  jim
 * Eliminated warnings from new generation of picky compilers.
 *
 * Revision 1.1010  1999/05/11 23:56:16  brunner
 * Changes for new charm version
 *
 * Revision 1.1009  1999/03/17 17:59:22  jim
 * Eliminated compiler warnings and errors.
 *
 * Revision 1.1008  1998/08/11 16:30:26  jim
 * Modified output from periodic boundary simulations to return atoms to
 * internally consistent coordinates.  We store the transformations which
 * were performed and undo them at the end.  It might be better to do this
 * by always keeping the original coordinates and only doing the transform
 * for the nonbonded terms but this works for now.
 *
 * Revision 1.1007  1998/03/03 23:05:02  brunner
 * Changed include files for new simplified Charm++ include file structure.
 *
 * Revision 1.1006  1998/01/15 05:40:52  jim
 * Added int return type to comparison operators.
 *
 * Revision 1.1005  1997/11/07 20:17:35  milind
 * Made NAMD to run on shared memory machines.
 *
 * Revision 1.1004  1997/07/09 21:26:39  milind
 * Ported NAMD2 to SP3. The SP specific code is within #ifdef SP2
 * and #endif's.
 *
 * Revision 1.1003  1997/03/19 11:53:59  ari
 * Add Broadcast mechanism.
 * Fixed RCS Log entries on files that did not have Log entries.
 * Added some register variables to Molecule and ComputeNonbondedExcl.C
 *
 ***************************************************************************/
