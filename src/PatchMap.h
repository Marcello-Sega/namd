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

class Patch;
class PatchMgr;
class Lattice;
class HomePatch;

class PatchMap
{
public:
  static PatchMap *Instance();
  inline static PatchMap *Object() { return _instance; }


  ~PatchMap(void);

  enum { MaxTwoAway = 5*5*5 - 3*3*3 };
  enum { MaxOneAway = 3*3*3 - 1 };
  enum { MaxOneOrTwoAway = MaxOneAway + MaxTwoAway };
  enum ErrCode { OK = 0, ERROR = -1 };

  static void registerPatchMgr(PatchMgr *pmgr) {
    patchMgr = pmgr;
  }

  HomePatchList *homePatchList();
  int numHomePatches(void);

  // numPatches() returns the number of patches being managed 
  // by the map.
  int numPatches(void);

  void setGridOriginAndLength(Vector o, Vector l);

  // assignToPatch(p,l) returns the patchID in which the position
  // p should be placed, based on the lattice l
  PatchID assignToPatch(Position p /* , Lattice &l */);

  // xDimension() returns the how many patches span the simulation space
  // along the x axis.
  int xDimension(void);

  // yDim() returns the how many patches span the simulation space
  // along the y axis.
  int yDimension(void);

  // zDim() returns the how many patches span the simulation space
  // along the z axis.
  int zDimension(void);

  int xIsPeriodic(void);
  int yIsPeriodic(void);
  int zIsPeriodic(void);

  // pid(xindex, yindex, zindex) returns the patch id for the given
  // patch coordinates.
  int pid(int xindex, int yindex, int zindex);

  // xIndex(pid) returns the x index for the given patch id.
  int xIndex(int pid);

  // yIndex(pid) returns the y index for the given patch id.
  int yIndex(int pid);

  // zIndex(pid) returns the z index for the given patch id.
  int zIndex(int pid);

  // node(pid) returns the node where the patch currently exists.
  int node(int pid);

  // minX(pid) returns the minimum x coordinate of the region of
  // space the patch is responsible for.
  Coordinate minX(int pid);

  // maxX(pid) returns the maximum x coordinate of the region of
  // space the patch is responsible for.
  Coordinate maxX(int pid);

  // minY(pid) returns the minimum y coordinate of the region of
  // space the patch is responsible for.
  Coordinate minY(int pid);

  // maxY(pid) returns the maximum y coordinate of the region of
  // space the patch is responsible for.
  Coordinate maxY(int pid);
  
  // minZ(pid) returns the minimum z coordinate of the region of
  // space the patch is responsible for.
  Coordinate minZ(int pid);

  // maxZ(pid) returns the maximum z coordinate of the region of
  // space the patch is responsible for.
  Coordinate maxZ(int pid);

  // numCids(pid) returns the number of compute ids which are registered
  // with the patch.  
  int numCids(int pid);
  
  // cid(pid,i) returns the i-th compute id registered
  // with the patch.  
  int cid(int pid, int i);

  void setPeriodicity(int x_per, int y_per, int z_per);

  // allocatePids(x_dim, y_dim, z_dim) tells the PatchMap to set aside
  // room for x_dim * y_dim * z_dim patches.
  ErrCode allocatePids(int x_dim, int y_dim, int z_dim);

  // PatchID requestPid(int *xi, int *yi, int *zi) tells the
  // PatchMap to return the next unused pid, and the indices of
  // that patch.
  PatchID requestPid(int *xi, int *yi, int *zi);

  // storePatch(pid, node, max_computes, x0, y0, z0, x1, y1, z1)
  // stores the info about the patch into the previously requested
  // pid.
  void storePatch(PatchID pid, int node, int max_computes,
		  Coordinate x0, Coordinate y0, Coordinate z0,
		  Coordinate x1, Coordinate y1, Coordinate z1);

  // newCid(pid,cid) stores a compute id associated with
  // patch id pid.  Error returned when there is no room to store
  // the pid.
  ErrCode newCid(int pid, int cid);

  // oneAwayNeighbors(pid, &n, neighbor_ids) returns the number 
  // and ids of adjacent patches.  The caller is expected to provide
  // sufficient storage for the neighbors.

  int oneAwayNeighbors(int pid, PatchID *neighbor_ids=0, int *transform_ids=0);

  // twoAwayNeighbors(pid, &n, neighbor_ids) returns the number 
  // and ids of all patches exactely two steps distant.  
  // The caller is expected to provide sufficient storage for the neighbors.

  int twoAwayNeighbors(int pid, PatchID *neighbor_ids, int *transform_ids = 0);

  int oneOrTwoAwayNeighbors(int pid, PatchID *neighbor_ids, int *transform_ids = 0);

  void printPatchMap(void);

  inline Patch *patch(PatchID pid);
  inline HomePatch *homePatch(PatchID pid);

  void registerPatch(PatchID pid, Patch *pptr);
  void unregisterPatch(PatchID pid, Patch *pptr);


protected:
  friend MapDistribMsg;
  void * pack (int *length);
  void unpack (void *in);
  
  PatchMap(void);
  
private:
  static PatchMap *_instance;
  static PatchMgr *patchMgr;
  
  struct PatchData
  {
    int node;
    int xi, yi, zi;
    Coordinate x0, x1, y0, y1, z0, z1;
    int numCids;
    int numCidsAllocated;
    ComputeID *cids;
    Patch *myPatch;
  };
  int curPatch;
  int nPatches;
  PatchData *patchData;
  int xDim, yDim, zDim;
  int xPeriodic, yPeriodic, zPeriodic;
  BigReal xOrigin, yOrigin, zOrigin;
  BigReal xLength, yLength, zLength;

};


//----------------------------------------------------------------------

inline Patch *PatchMap::patch(PatchID pid)
{
  return patchData[pid].myPatch;
}

inline HomePatch *PatchMap::homePatch(PatchID pid)
{
  return (HomePatch *)patchData[pid].myPatch;
}

#endif /* PATCHMAP_H */


/***************************************************************************
 * RCS INFORMATION:
 *
 *	$RCSfile: PatchMap.h,v $
 *	$Author: ari $	$Locker:  $		$State: Exp $
 *	$Revision: 1.1000 $	$Date: 1997/02/06 15:59:04 $
 *
 ***************************************************************************
 * REVISION HISTORY:
 *
 * $Log: PatchMap.h,v $
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
