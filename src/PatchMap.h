//-*-c++-*-
/***************************************************************************/
/*                                                                         */
/*              (C) Copyright 1996 The Board of Trustees of the            */
/*                          University of Illinois                         */
/*                           All Rights Reserved                           */
/*									   */
/***************************************************************************/

/***************************************************************************
 * DESCRIPTION: PatchMap.h
 *
 ***************************************************************************/

#ifndef PATCHMAP_H
#define PATCHMAP_H

#include "NamdTypes.h"
#include "HomePatchList.h"
#include "Lattice.h"
#include "ProcessorPrivate.h"

class Patch;
class PatchMgr;
class HomePatch;

class PatchMap
{
public:
  static PatchMap *Instance();
  inline static PatchMap *Object() { return CpvAccess(PatchMap_instance); }

  void initialize(ScaledPosition xmin, ScaledPosition xmax,
                      Lattice lattice, BigReal patchSize);
  void checkMap();

  ~PatchMap(void);

  enum { MaxTwoAway = 5*5*5 - 3*3*3 };
  enum { MaxOneAway = 3*3*3 - 1 };
  enum { MaxOneOrTwoAway = MaxOneAway + MaxTwoAway };

  static void registerPatchMgr(PatchMgr *pmgr) {
    CpvAccess(PatchMap_patchMgr) = pmgr;
  }

  HomePatchList *homePatchList();
  int numHomePatches(void);

  // returns the number of patches being managed 
  inline int numPatches(void) const { return nPatches; }

  // returns the number of patches in each dimension
  inline int gridsize_a(void) const { return aDim; }
  inline int gridsize_b(void) const { return bDim; }
  inline int gridsize_c(void) const { return cDim; }

  // returns 1 if periodic in each dimension
  inline int periodic_a(void) const { return aPeriodic; }
  inline int periodic_b(void) const { return bPeriodic; }
  inline int periodic_c(void) const { return cPeriodic; }

  // returns the origin (minimum, not center) of patch grid
  inline ScaledPosition origin(void) const {
    return ScaledPosition(aOrigin,bOrigin,cOrigin);
  }

  // returns the patch id for the given indices
  inline int pid(int aIndex, int bIndex, int cIndex);

  // returns the [abc] index for the given patch id.
  inline int index_a(int pid) const { return pid % aDim; }
  inline int index_b(int pid) const { return (pid / aDim) % bDim; }
  inline int index_c(int pid) const { return pid / (aDim*bDim); }

  // returns the min/max [abc] scaled coordinate
  inline BigReal min_a(int pid) const { return patchData[pid].aMin; }
  inline BigReal max_a(int pid) const { return patchData[pid].aMax; }
  inline BigReal min_b(int pid) const { return patchData[pid].bMin; }
  inline BigReal max_b(int pid) const { return patchData[pid].bMax; }
  inline BigReal min_c(int pid) const { return patchData[pid].cMin; }
  inline BigReal max_c(int pid) const { return patchData[pid].cMax; }

  // asssigns atom to patch based on position and lattice
  inline PatchID assignToPatch(Position p, const Lattice &l);

  // gives more downstream patch of pid1, pid2; handles periodicity right
  // given patches must be neighbors!!!
  inline int downstream(int pid1, int pid2);

  // returns the node where the patch currently exists.
  inline int node(int pid) const { return patchData[pid].node; }

  // numCids(pid) returns the number of compute ids which are registered
  inline int numCids(int pid) const { return patchData[pid].numCids; }
  
  // cid(pid,i) returns the i-th compute id registered
  inline int cid(int pid, int i) const { return patchData[pid].cids[i]; }

  void assignNode(PatchID, NodeID);

  // newCid(pid,cid) stores a compute id associated with
  // patch id pid.  Error returned when there is no room to store
  // the pid.
  void newCid(int pid, int cid);

  // oneAwayNeighbors(pid, &n, neighbor_ids) returns the number 
  // and ids of adjacent patches.  The caller is expected to provide
  // sufficient storage for the neighbors.

  int oneAwayNeighbors(int pid, PatchID *neighbor_ids=0, int *transform_ids=0);

  // twoAwayNeighbors(pid, &n, neighbor_ids) returns the number 
  // and ids of all patches exactely two steps distant.  
  // The caller is expected to provide sufficient storage for the neighbors.

  int twoAwayNeighbors(int pid, PatchID *neighbor_ids, int *transform_ids = 0);

  int oneOrTwoAwayNeighbors(int pid, PatchID *neighbor_ids,
			    int *transform_ids = 0);

  int upstreamNeighbors(int pid, PatchID *neighbor_ids, 
			int *transform_ids = 0);

  int downstreamNeighbors(int pid, PatchID *neighbor_ids, 
			  int *transform_ids = 0);

  void printPatchMap(void);

  inline Patch *patch(PatchID pid);
  HomePatch *homePatch(PatchID pid);

  void registerPatch(PatchID pid, HomePatch *pptr);
  void unregisterPatch(PatchID pid, HomePatch *pptr);

  void registerPatch(PatchID pid, Patch *pptr);
  void unregisterPatch(PatchID pid, Patch *pptr);


protected:
  friend class MapDistribMsg;
  int packSize(void);
  void pack(char *buf);
  void unpack(char *buf);
  
  PatchMap(void);
  
private:
  struct PatchData
  {
    int node;
    int aIndex, bIndex, cIndex;
    Coordinate aMin, aMax, bMin, bMax, cMin, cMax;
    int numCids;
    int numCidsAllocated;
    ComputeID *cids;
    Patch *myPatch;
    HomePatch *myHomePatch;
  };
  int nPatches;
  PatchData *patchData;
  int aDim, bDim, cDim;
  int aPeriodic, bPeriodic, cPeriodic;
  int aMaxIndex, bMaxIndex, cMaxIndex;
  BigReal aOrigin, bOrigin, cOrigin;
  BigReal aLength, bLength, cLength;

};


//----------------------------------------------------------------------

inline Patch *PatchMap::patch(PatchID pid)
{
  return patchData[pid].myPatch;
}

#endif /* PATCHMAP_H */


/***************************************************************************
 * RCS INFORMATION:
 *
 *	$RCSfile: PatchMap.h,v $
 *	$Author: jim $	$Locker:  $		$State: Exp $
 *	$Revision: 1.1014 $	$Date: 1999/09/03 20:46:19 $
 *
 ***************************************************************************
 * REVISION HISTORY:
 *
 * $Log: PatchMap.h,v $
 * Revision 1.1014  1999/09/03 20:46:19  jim
 * Support for non-orthogonal periodic boundary conditions.
 *
 * Revision 1.1013  1999/08/11 16:52:23  jim
 * Make homePatch() method return NULL for proxies rather than casting.
 *
 * Revision 1.1012  1999/05/11 23:56:42  brunner
 * Changes for new charm version
 *
 * Revision 1.1011  1998/07/16 18:52:13  jim
 * Localized common downstream patch optimization.
 *
 * Revision 1.1010  1998/06/24 23:40:33  brunner
 * Added downstreamNeighbors() and LdbCoordinator fixes.  I don't know
 * why Patch.C is different.
 *
 * Revision 1.1009  1998/06/18 14:44:36  jim
 * Eliminated warnings and errors from aCC.
 *
 * Revision 1.1008  1997/11/07 20:17:44  milind
 * Made NAMD to run on shared memory machines.
 *
 * Revision 1.1007  1997/10/06 00:12:33  jim
 * Added PatchMap.inl, sped up cycle-boundary tuple code.
 *
 * Revision 1.1006  1997/09/28 22:36:53  jim
 * Modified tuple-based computations to not duplicate calculations and
 * only require "upstream" proxies.
 *
 * Revision 1.1005  1997/04/10 09:14:02  ari
 * Final debugging for compute migration / proxy creation for load balancing.
 * Lots of debug code added, mostly turned off now.
 * Fixed bug in PositionBox when Patch had no dependencies.
 * Eliminated use of cout and misuse of iout in numerous places.
 *                                            Ari & Jim
 *
 * Revision 1.1004  1997/03/27 08:04:21  jim
 * Reworked Lattice to keep center of cell fixed during rescaling.
 *
 * Revision 1.1003  1997/03/14 21:40:14  ari
 * Reorganized startup to make possible inital load
 * balancing by changing methods in WorkDistrib.
 * Also made startup more transparent and easier
 * to modify.
 *
 * Revision 1.1002  1997/02/13 04:43:12  jim
 * Fixed initial hanging (bug in PatchMap, but it still shouldn't have
 * happened) and saved migration messages in the buffer from being
 * deleted, but migration still dies (even on one node).
 *
 * Revision 1.1001  1997/02/07 16:56:52  nealk
 * Added Origin() to return origin.
 *
 * Revision 1.1000  1997/02/06 15:59:04  ari
 * Resetting CVS to merge branches back into the main trunk.
 * We will stick to main trunk development as suggested by CVS manual.
 * We will set up tags to track fixed points of development/release
 * as suggested by CVS manual - all praise the CVS manual.
 *
 * Revision 1.779  1997/02/06 15:53:22  ari
 * Updating Revision Line, getting rid of branches
 *
 * Revision 1.778.2.3  1997/02/06 02:35:31  jim
 * Implemented periodic boundary conditions - may not work with
 * atom migration yet, but doesn't seem to alter calculation,
 * appears to work correctly when turned on.
 * NamdState chdir's to same directory as config file in argument.
 *
 * Revision 1.778.2.2  1997/02/05 22:18:19  ari
 * Added migration code - Currently the framework is
 * there with compiling code.  This version does
 * crash shortly after migration is complete.
 * Migration appears to complete, but Patches do
 * not appear to be left in a correct state.
 *
 * Revision 1.778.2.1  1997/01/28 17:28:50  jim
 * First top-down changes for periodic boundary conditions, added now to
 * avoid conflicts with Ari's migration system.
 *
 * Revision 1.778  1997/01/28 00:31:11  ari
 * internal release uplevel to 1.778
 *
 * Revision 1.777.2.2  1997/01/27 22:45:33  ari
 * Basic Atom Migration Code added.
 * Added correct magic first line to .h files for xemacs to go to C++ mode.
 * Compiles and runs without migration turned on.
 *
 * Revision 1.777.2.1  1997/01/21 23:04:48  ari
 * Basic framework for atom migration placed into code.  - Non
 * functional since it is not called.  Works currently without
 * atom migration.
 *
 * Revision 1.777  1997/01/17 19:36:47  ari
 * Internal CVS leveling release.  Start development code work
 * at 1.777.1.1.
 *
 * Revision 1.10  1997/01/13 21:04:19  jim
 * added numHomePatches()
 *
 * Revision 1.9  1996/12/18 21:07:54  jim
 * added oneOrTwoAwayNeighbors()
 *
 * Revision 1.8  1996/12/12 08:57:17  jim
 * added MapDistribMsg packing / unpacking routines
 *
 * Revision 1.7  1996/11/21 20:39:29  jim
 * small bug fixes
 *
 * Revision 1.6  1996/11/01 21:20:45  ari
 * *** empty log message ***
 *
 * Revision 1.5  1996/10/29 23:35:27  ari
 * *** empty log message ***
 *
 * Revision 1.4  1996/10/10 17:23:24  brunner
 * Added patch * in patchmap
 *
 * Revision 1.3  1996/08/23 22:03:52  brunner
 * *** empty log message ***
 *
 * Revision 1.2  1996/08/19 21:37:02  brunner
 * Changed Position to Coordinate
 *
 * Revision 1.1  1996/08/16 20:43:53  brunner
 * Initial revision
 *
 * Revision 1.7  1996/08/03 20:08:09  brunner
 * *** empty log message ***
 *
 * Revision 1.6  1996/07/16 20:06:24  brunner
 * *** empty log message ***
 *
 * Revision 1.5  1996/07/16 19:59:00  brunner
 * *** empty log message ***
 *
 * Revision 1.2  1996/06/11 22:36:35  brunner
 * *** empty log message ***
 *
 * Revision 1.1  1996/06/10 18:52:26  brunner
 * Initial revision
 *
 *
 ***************************************************************************/
