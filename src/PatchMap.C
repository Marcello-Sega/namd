/***************************************************************************/
/*         (C) Copyright 1996,1997 The Board of Trustees of the            */
/*                          University of Illinois                         */
/*                           All Rights Reserved                           */
/***************************************************************************/
/***************************************************************************
 * DESCRIPTION:
 *
 ***************************************************************************/

#include <stddef.h>
#include <unistd.h>
#include <stdio.h>

#include "chare.h"
#include "ckdefs.h"
#include "c++interface.h"

#include "PatchMgr.h"
#include "PatchMap.h"
#include "Patch.h"
#include "Lattice.h"
#include "HomePatchList.h"

#define DEBUGM
#define MIN_DEBUG_LEVEL 5
#include "Debug.h"

// static initialization
PatchMap *PatchMap::_instance = 0;

// Safe singleton creation
PatchMap *PatchMap::Instance() {
  if (_instance == 0) {
     _instance = new PatchMap;
  }
  return(_instance);
}


PatchMap::PatchMap(void)
{
  nPatches = 0;
  patchData = NULL;
  xDim = yDim = zDim = 0;
  xPeriodic = yPeriodic = zPeriodic = 0;
}

PatchMap::~PatchMap(void)
{
  if (patchData)
  {
    int i;

    for (i=0; i<nPatches;i++)
    {
      delete [] patchData[i].cids;
      patchData[i].cids=NULL;
    }
    delete [] patchData;
    patchData=NULL;
    nPatches=0;
  }
  
}

#undef PACK
#define PACK(type,data) { *((type*)b) = data; b += sizeof(type); }

void * PatchMap::pack (int *length)
{
  DebugM(4,"Packing PatchMap on node " << CMyPe() << endl);
  int i,j;

  // calculate memory needed
  int size = 0;
  size += 8 * sizeof(int) + 6 * sizeof(BigReal);
  for(i=0;i<nPatches;++i)
  {
    size += sizeof(PatchData);
    size += patchData[i].numCidsAllocated * sizeof(ComputeID);
  }
  *length = size;

  // allocate needed memory
  char * const buffer = new char[size];

  // fill in the data
  char *b = buffer;
  PACK(int,curPatch);
  PACK(int,nPatches);
  DebugM(4,"nPatches = " << nPatches << endl);
  PACK(int,xDim); PACK(int,yDim); PACK(int,zDim);
  PACK(int,xPeriodic); PACK(int,yPeriodic); PACK(int,zPeriodic);
  PACK(BigReal,xOrigin); PACK(BigReal,yOrigin); PACK(BigReal,zOrigin);
  PACK(BigReal,xLength); PACK(BigReal,yLength); PACK(BigReal,zLength);
  for(i=0;i<nPatches;++i)
  {
    DebugM(3,"Packing Patch " << i << " is on node " << patchData[i].node << 
	" with " << patchData[i].numCidsAllocated << " allocated.\n");
    PACK(PatchData,patchData[i]);
    for(j=0;j<patchData[i].numCidsAllocated;++j)
      PACK(ComputeID,patchData[i].cids[j]);
  }
  DebugM(3,buffer + size - b << " == 0 ?" << endl);

  return buffer;
}

#undef UNPACK
#define UNPACK(type,data) { data = *((type*)b); b += sizeof(type); }

void PatchMap::unpack (void *in)
{
  DebugM(4,"Unpacking PatchMap on node " << CMyPe() << endl);
  int i,j;
  char *b = (char*)in;
  UNPACK(int,curPatch);
  UNPACK(int,nPatches);
  DebugM(4,"nPatches = " << nPatches << endl);
  UNPACK(int,xDim); UNPACK(int,yDim); UNPACK(int,zDim);
  UNPACK(int,xPeriodic); UNPACK(int,yPeriodic); UNPACK(int,zPeriodic);
  UNPACK(BigReal,xOrigin); UNPACK(BigReal,yOrigin); UNPACK(BigReal,zOrigin);
  UNPACK(BigReal,xLength); UNPACK(BigReal,yLength); UNPACK(BigReal,zLength);
  patchData = new PatchData[nPatches];
  for(i=0;i<nPatches;++i)
  {
    UNPACK(PatchData,patchData[i]);
    DebugM(3,"Unpacking Patch " << i << " is on node " << patchData[i].node << 
	" with " << patchData[i].numCidsAllocated << " allocated.\n");
    patchData[i].cids = new ComputeID[patchData[i].numCidsAllocated];
    for(j=0;j<patchData[i].numCidsAllocated;++j)
      UNPACK(ComputeID,patchData[i].cids[j]);
  }
}


//---------------------------------------------------------------------
// Access HomePatch information

int PatchMap::numHomePatches(void)
{
  return patchMgr->homePatches.size();
}

HomePatchList *PatchMap::homePatchList() {
  return &(patchMgr->homePatches);
}

//----------------------------------------------------------------------
int PatchMap::numPatches(void)
{
  return nPatches;
}

//----------------------------------------------------------------------
void PatchMap::setGridOriginAndLength(Vector o, Vector l)
{
  xOrigin = o.x; yOrigin = o.y; zOrigin = o.z;
  xLength = l.x; yLength = l.y; zLength = l.z;
}

//----------------------------------------------------------------------
PatchID PatchMap::assignToPatch(Position p)
{
  int xi, yi, zi;
  xi = (int)floor(((BigReal)xDim)*((p.x-xOrigin)/xLength));
  yi = (int)floor(((BigReal)yDim)*((p.y-yOrigin)/yLength));
  zi = (int)floor(((BigReal)zDim)*((p.z-zOrigin)/zLength));
  return pid(xi,yi,zi);
}

//----------------------------------------------------------------------
int PatchMap::xDimension(void)
{
  return xDim;
}

//----------------------------------------------------------------------
int PatchMap::yDimension(void)
{
  return yDim;
}

//----------------------------------------------------------------------
int PatchMap::zDimension(void)
{
  return zDim;
}

//----------------------------------------------------------------------
Vector PatchMap::Origin(void)
{
  Vector o;
  o.x = xOrigin;
  o.y = yOrigin;
  o.z = zOrigin;
  return(o);
}

//----------------------------------------------------------------------
int PatchMap::xIsPeriodic(void)
{
  return xPeriodic;
}

//----------------------------------------------------------------------
int PatchMap::yIsPeriodic(void)
{
  return yPeriodic;
}

//----------------------------------------------------------------------
int PatchMap::zIsPeriodic(void)
{
  return zPeriodic;
}

//----------------------------------------------------------------------
#define MODULO(I,J) ( (I)<0 ? (I)-(J)*((I)/(J)-1) : (I)-(J)*((I)/(J)) )

int PatchMap::pid(int xIndex, int yIndex, int zIndex)
{
  int allsame = 0;
  if ( xPeriodic ) xIndex = MODULO(xIndex,xDim);
  else
  {
    if ( xIndex < 0 ) xIndex = 0;
    if ( xIndex >= xDim ) xIndex = xDim - 1;
  }
  if ( yPeriodic ) yIndex = MODULO(yIndex,yDim);
  else
  {
    if ( yIndex < 0 ) yIndex = 0;
    if ( yIndex >= yDim ) yIndex = yDim - 1;
  }
  if ( zPeriodic ) zIndex = MODULO(zIndex,zDim);
  else
  {
    if ( zIndex < 0 ) zIndex = 0;
    if ( zIndex >= zDim ) zIndex = zDim - 1;
  }
  return ((zIndex*yDim)+yIndex)*xDim + xIndex;
}

//----------------------------------------------------------------------
int PatchMap::xIndex(int pid)
{
  return pid % xDim;
}

//----------------------------------------------------------------------
int PatchMap::yIndex(int pid)
{
  return (pid / xDim) % yDim;
}

//----------------------------------------------------------------------
int PatchMap::zIndex(int pid)
{
  return (pid / (xDim*yDim));
}

//----------------------------------------------------------------------
int PatchMap::node(int pid)
{
  return patchData[pid].node;
}

//----------------------------------------------------------------------
Coordinate PatchMap::minX(int pid)
{
  return patchData[pid].x0;
}

//----------------------------------------------------------------------
Coordinate PatchMap::maxX(int pid)
{
  return patchData[pid].x1;
}

//----------------------------------------------------------------------
Coordinate PatchMap::minY(int pid)
{
  return patchData[pid].y0;
}

//----------------------------------------------------------------------
Coordinate PatchMap::maxY(int pid)
{
  return patchData[pid].y1;
}

//----------------------------------------------------------------------
Coordinate PatchMap::minZ(int pid)
{
  return patchData[pid].z0;
}

//----------------------------------------------------------------------
Coordinate PatchMap::maxZ(int pid)
{
  return patchData[pid].z1;
}

//----------------------------------------------------------------------
int PatchMap::numCids(int pid)
{
  return patchData[pid].numCids;
}

//----------------------------------------------------------------------
int PatchMap::cid(int pid,int i)
{
  return patchData[pid].cids[i];
}

//----------------------------------------------------------------------
void PatchMap::setPeriodicity(int x_per, int y_per, int z_per)
{
  xPeriodic = x_per;
  yPeriodic = y_per;
  zPeriodic = z_per;
}

//----------------------------------------------------------------------
PatchMap::ErrCode PatchMap::allocatePids(int ixDim, int iyDim, int izDim)
{
  int i;

  if (patchData)
  {
    for (i=0; i<nPatches;i++)
    {
      delete [] patchData[i].cids;
      patchData[i].cids=NULL;
    }
    delete [] patchData;
  }
  curPatch=0;
  xDim = ixDim;
  yDim = iyDim;
  zDim = izDim;
  nPatches=xDim*yDim*zDim;
  patchData = new PatchData[nPatches];
  if (!patchData)
    return ERROR;

  for(i=0;i<nPatches;i++)
  {
    patchData[i].numCids=0;
    patchData[i].xi = xIndex(i);
    patchData[i].yi = yIndex(i);
    patchData[i].zi = zIndex(i);
    patchData[i].cids=NULL;
    patchData[i].numCidsAllocated=0;
    patchData[i].myPatch=NULL;
  }

  return OK;
}

//----------------------------------------------------------------------
PatchID PatchMap::requestPid(int *xi, int *yi, int *zi)
{
  int pid;
  
  if ((patchData != NULL) && (curPatch < nPatches))
  {
    pid = curPatch;
    curPatch++;
    *xi = patchData[pid].xi;
    *yi = patchData[pid].yi;
    *zi = patchData[pid].zi;
    return pid;
  }
  return -1;
}

//----------------------------------------------------------------------
void PatchMap::storePatchCoord(PatchID pid,
				ScaledPosition min, ScaledPosition max)
{
  patchData[pid].x0 = min.x;
  patchData[pid].x1 = max.x;
  patchData[pid].y0 = min.y;
  patchData[pid].y1 = max.y;
  patchData[pid].z0 = min.z;
  patchData[pid].z1 = max.z;
}

void PatchMap::assignNode(PatchID pid, NodeID node) {
  patchData[pid].node=node;
}

void PatchMap::allocateCompute(PatchID pid, int max_computes) {
  patchData[pid].numCids = 0;
  patchData[pid].cids = new int[max_computes];
  for ( int i = 0; i < max_computes; ++i ) patchData[pid].cids[i] = -1;
  patchData[pid].numCidsAllocated = max_computes;
}

//----------------------------------------------------------------------
PatchMap::ErrCode PatchMap::newCid(int pid, int cid)
{
  if (patchData[pid].numCids < patchData[pid].numCidsAllocated)
  {
    patchData[pid].cids[patchData[pid].numCids]=cid;
    patchData[pid].numCids++;
    return OK;
  } else return ERROR;
}

//----------------------------------------------------------------------
int PatchMap::oneAwayNeighbors(int pid, PatchID *neighbor_ids, int *transform_ids)
{
  int xi, yi, zi;
  int xinc, yinc, zinc;
  int n=0;

  for(zinc=-1;zinc<=1;zinc++)
  {
    zi = patchData[pid].zi + zinc;
    if ((zi < 0) || (zi >= zDim))
      if ( ! zPeriodic ) continue;
    for(yinc=-1;yinc<=1;yinc++)
    {
      yi = patchData[pid].yi + yinc;
      if ((yi < 0) || (yi >= yDim))
	if ( ! yPeriodic ) continue;
      for(xinc=-1;xinc<=1;xinc++)
      {
	if ((xinc==0) && (yinc==0) && (zinc==0))
	  continue;

	xi = patchData[pid].xi + xinc;
	if ((xi < 0) || (xi >= xDim))
	  if ( ! xPeriodic ) continue;

	if (neighbor_ids)
	  neighbor_ids[n]=this->pid(xi,yi,zi);
	if ( transform_ids )
	{
	  int xt = 0; if ( xi < 0 ) xt = -1; if ( xi >= xDim ) xt = 1;
	  int yt = 0; if ( yi < 0 ) yt = -1; if ( yi >= yDim ) yt = 1;
	  int zt = 0; if ( zi < 0 ) zt = -1; if ( zi >= zDim ) zt = 1;
	  transform_ids[n] = Lattice::index(xt,yt,zt);
	}
	n++;
      }
    }
  }
  DebugM(3,"Patch " << pid << " has " << n << " first neighbors.\n");
  return n;
}

//----------------------------------------------------------------------
int PatchMap::twoAwayNeighbors(int pid, PatchID *neighbor_ids,  int *transform_ids)
{
  // ELIMINATE USE OF TWO-AWAY NEIGHBORS

  return 0;

  // ELIMINATE USE OF TWO-AWAY NEIGHBORS

  /*
  int xi, yi, zi;
  int xinc, yinc, zinc;
  int n=0;

  for(zinc=-2;zinc<=2;zinc++)
  {
    zi = patchData[pid].zi + zinc;
    if ((zi < 0) || (zi >= zDim))
      if ( ! zPeriodic ) continue;
    for(yinc=-2;yinc<=2;yinc++)
    {
      yi = patchData[pid].yi + yinc;
      if ((yi < 0) || (yi >= yDim))
	if ( ! yPeriodic ) continue;
      for(xinc=-2;xinc<=2;xinc++)
      {
	if (!((xinc==2) || (yinc==2) || (zinc==2) ||
	    (xinc==-2) || (yinc==-2) || (zinc==-2)))
	  continue;

	xi = patchData[pid].xi + xinc;
	if ((xi < 0) || (xi >= xDim))
	  if ( ! xPeriodic ) continue;

	neighbor_ids[n]=this->pid(xi,yi,zi);
	if ( transform_ids )
	{
	  int xt = 0; if ( xi < 0 ) xt = 1; if ( xi >= xDim ) xt = -1;
	  int yt = 0; if ( yi < 0 ) yt = 1; if ( yi >= yDim ) yt = -1;
	  int zt = 0; if ( zi < 0 ) zt = 1; if ( zi >= zDim ) zt = -1;
	  transform_ids[n] = Lattice::index(xt,yt,zt);
	}
	n++;
      }
    }
  }
  DebugM(3,"Patch " << pid << " has " << n << " second neighbors.\n");
  return n;
  */
}

//----------------------------------------------------------------------
int PatchMap::oneOrTwoAwayNeighbors(int pid, PatchID *neighbor_ids, int *transform_ids)
{
  int numOneAway = oneAwayNeighbors(pid,neighbor_ids,transform_ids);
  int numTwoAway = twoAwayNeighbors(pid,neighbor_ids+numOneAway,
			transform_ids?transform_ids+numOneAway:0);
  return numOneAway + numTwoAway;
}

//----------------------------------------------------------------------
void PatchMap::printPatchMap(void)
{
  CPrintf("---------------------------------------");
  CPrintf("---------------------------------------\n");

  CPrintf("curPatch = %d\n",curPatch);  
  CPrintf("nPatches = %d\n",nPatches);
  for(int i=0;i<nPatches;i++)
  {
    CPrintf("Patch %d:\n",i);
    CPrintf("  node = %d\n",patchData[i].node);
    CPrintf("  xi,yi,zi = %d, %d, %d\n",
	    patchData[i].xi,patchData[i].yi,patchData[i].zi);
    CPrintf("  x0,y0,z0 = %f, %f, %f\n",
	    patchData[i].x0,patchData[i].y0,patchData[i].z0);
    CPrintf("  x1,y1,z1 = %f, %f, %f\n",
	    patchData[i].x1,patchData[i].y1,patchData[i].z1);
    CPrintf("  numCids = %d\n",patchData[i].numCids);
    CPrintf("  numCidsAllocated = %d\n",patchData[i].numCidsAllocated);
    for(int j=0; j < patchData[i].numCids; j++)
    {
      CPrintf(" %10d ",patchData[i].cids[j]);
      if (!((j+1) % 6))
	CPrintf("\n");
    }
    CPrintf("\n---------------------------------------");
    CPrintf("---------------------------------------\n");
  }

}

//----------------------------------------------------------------------
void PatchMap::registerPatch(PatchID pid, Patch *pptr)
{
  patchData[pid].myPatch = pptr;
}

//----------------------------------------------------------------------
void PatchMap::unregisterPatch(PatchID pid, Patch *pptr)
{
  if (pptr == patchData[pid].myPatch)
      patchData[pid].myPatch = NULL;
}



/***************************************************************************
 * RCS INFORMATION:
 *
 *	$RCSfile: PatchMap.C,v $
 *	$Author: jim $	$Locker:  $		$State: Exp $
 *	$Revision: 1.1008 $	$Date: 1997/03/27 08:04:20 $
 *
 ***************************************************************************
 * REVISION HISTORY:
 *
 * $Log: PatchMap.C,v $
 * Revision 1.1008  1997/03/27 08:04:20  jim
 * Reworked Lattice to keep center of cell fixed during rescaling.
 *
 * Revision 1.1007  1997/03/19 11:54:47  ari
 * Add Broadcast mechanism.
 * Fixed RCS Log entries on files that did not have Log entries.
 * Added some register variables to Molecule and ComputeNonbondedExcl.C
 *
 * Revision 1.1006  1997/03/14 21:40:14  ari
 * Reorganized startup to make possible inital load
 * balancing by changing methods in WorkDistrib.
 * Also made startup more transparent and easier
 * to modify.
 *
 * Revision 1.1005  1997/02/28 23:14:23  jim
 * Eliminated use of two-away neighbors, method now returns 0.
 *
 * Revision 1.1004  1997/02/13 04:43:11  jim
 * Fixed initial hanging (bug in PatchMap, but it still shouldn't have
 * happened) and saved migration messages in the buffer from being
 * deleted, but migration still dies (even on one node).
 *
 * Revision 1.1003  1997/02/07 22:52:16  jim
 * Eliminated use of nAtomBased and uninitialized memory reads.
 *
 * Revision 1.1002  1997/02/07 17:39:40  ari
 * More debugging for atomMigration.
 * Using -w on CC got us some minor fixes
 * using purify got us a major memory problem due to bad sizing of dummy force
 *
 * Revision 1.1001  1997/02/07 16:56:50  nealk
 * Added Origin() to return origin.
 *
 * Revision 1.1000  1997/02/06 15:59:03  ari
 * Resetting CVS to merge branches back into the main trunk.
 * We will stick to main trunk development as suggested by CVS manual.
 * We will set up tags to track fixed points of development/release
 * as suggested by CVS manual - all praise the CVS manual.
 *
 * Revision 1.779  1997/02/06 15:53:22  ari
 * Updating Revision Line, getting rid of branches
 *
 * Revision 1.778.2.3  1997/02/06 02:35:30  jim
 * Implemented periodic boundary conditions - may not work with
 * atom migration yet, but doesn't seem to alter calculation,
 * appears to work correctly when turned on.
 * NamdState chdir's to same directory as config file in argument.
 *
 * Revision 1.778.2.2  1997/02/05 22:18:18  ari
 * Added migration code - Currently the framework is
 * there with compiling code.  This version does
 * crash shortly after migration is complete.
 * Migration appears to complete, but Patches do
 * not appear to be left in a correct state.
 *
 * Revision 1.778.2.1  1997/01/28 17:28:49  jim
 * First top-down changes for periodic boundary conditions, added now to
 * avoid conflicts with Ari's migration system.
 *
 * Revision 1.778  1997/01/28 00:31:10  ari
 * internal release uplevel to 1.778
 *
 * Revision 1.777  1997/01/17 19:36:47  ari
 * Internal CVS leveling release.  Start development code work
 * at 1.777.1.1.
 *
 * Revision 1.15  1997/01/13 21:04:19  jim
 * added numHomePatches()
 *
 * Revision 1.14  1996/12/19 00:37:31  jim
 * increase MIN_DEBUG_LEVEL
 *
 * Revision 1.13  1996/12/18 21:07:54  jim
 * added oneOrTwoAwayNeighbors()
 *
 * Revision 1.12  1996/12/13 22:58:01  nealk
 * Found pack/unpack bug and corrected it.  (wrong offset!)
 *
 * Revision 1.11  1996/12/13 19:39:55  jim
 * added debugging, looking for error in PatchMap sending
 *
 * Revision 1.10  1996/12/12 08:57:17  jim
 * added MapDistribMsg packing / unpacking routines
 *
 * Revision 1.9  1996/11/21 21:17:10  jim
 * small bug fix in patch indexing functions
 *
 * Revision 1.8  1996/11/01 21:20:45  ari
 * *** empty log message ***
 *
 * Revision 1.7  1996/10/29 23:35:27  ari
 * *** empty log message ***
 *
 * Revision 1.6  1996/10/16 08:22:39  ari
 * *** empty log message ***
 *
 * Revision 1.5  1996/10/10 17:23:41  brunner
 * Added patch * in patchmap
 *
 * Revision 1.4  1996/08/23 22:03:52  brunner
 * *** empty log message ***
 *
 * Revision 1.3  1996/08/19 21:37:02  brunner
 * Changed Position to Coordinate
 *
 * Revision 1.2  1996/08/16 21:41:11  brunner
 * *** empty log message ***
 *
 * Revision 1.1  1996/08/16 20:43:53  brunner
 * Initial revision
 *
 * Revision 1.2  1996/08/03 20:08:09  brunner
 * *** empty log message ***
 *
 * Revision 1.1  1996/07/16 20:07:08  brunner
 * Initial revision
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
