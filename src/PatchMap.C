/***************************************************************************/
/*                                                                         */
/*              (C) Copyright 1996 The Board of Trustees of the            */
/*                          University of Illinois                         */
/*                           All Rights Reserved                           */
/*									   */
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

#include "PatchMap.h"
#include "Patch.h"

#define DEBUGM
#define MIN_DEBUG_LEVEL 3
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
  int xDim = yDim = zDim = 0;
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
  DebugM(4,"Packing PatchMap\n");
  int i,j;

  // calculate memory needed
  int size = 0;
  size += 5 * sizeof(int);
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
  PACK(int,xDim); PACK(int,yDim); PACK(int,zDim);
  for(i=0;i<nPatches;++i)
  {
    PACK(PatchData,patchData[i]);
    for(j=0;j<patchData[i].numCidsAllocated;++j)
      PACK(ComputeID,patchData[i].cids[j]);
  }

  return buffer;
}

#undef UNPACK
#define UNPACK(type,data) { data = *((type*)b); b += sizeof(type); }

void PatchMap::unpack (void *in)
{
  DebugM(4,"Unpacking PatchMap\n");
  int i,j;
  char *b = (char*)in;
  UNPACK(int,curPatch);
  UNPACK(int,nPatches);
  UNPACK(int,xDim); UNPACK(int,yDim); UNPACK(int,zDim);
  patchData = new PatchData[nPatches];
  for(i=0;i<nPatches;++i)
  {
    UNPACK(PatchData,patchData[i]);
    patchData[i].cids = new ComputeID[patchData[i].numCidsAllocated];
    for(j=0;j<patchData[i].numCidsAllocated;++j)
      UNPACK(ComputeID,patchData[i].cids[j]);
  }
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
int PatchMap::pid(int xIndex, int yIndex, int zIndex)
{
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
void PatchMap::storePatch(PatchID pid, int node, int max_computes,
			  Coordinate x0, Coordinate y0, Coordinate z0,
			  Coordinate x1, Coordinate y1, Coordinate z1)
{
  patchData[pid].node=node;
  patchData[pid].numCids = 0;
  patchData[pid].cids = new int[max_computes];
  patchData[pid].numCidsAllocated = max_computes;
  patchData[pid].x0 = x0;
  patchData[pid].x1 = x1;
  patchData[pid].y0 = y0;
  patchData[pid].y1 = y1;
  patchData[pid].z0 = z0;
  patchData[pid].z1 = z1;
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
int PatchMap::oneAwayNeighbors(int pid, PatchID *neighbor_ids)
{
  int xi, yi, zi;
  int xinc, yinc, zinc;
  int n=0;

  for(zinc=-1;zinc<=1;zinc++)
  {
    zi = patchData[pid].zi + zinc;
    if ((zi < 0) || (zi >= zDim))
      continue;
    for(yinc=-1;yinc<=1;yinc++)
    {
      yi = patchData[pid].yi + yinc;
      if ((yi < 0) || (yi >= yDim))
	continue;
      for(xinc=-1;xinc<=1;xinc++)
      {
	if ((xinc==0) && (yinc==0) && (zinc==0))
	  continue;

	xi = patchData[pid].xi + xinc;
	if ((xi < 0) || (xi >= xDim))
	  continue;

	neighbor_ids[n]=this->pid(xi,yi,zi);
	n++;
      }
    }
  }
  return n;
}

//----------------------------------------------------------------------
int PatchMap::twoAwayNeighbors(int pid, PatchID *neighbor_ids)
{
  int xi, yi, zi;
  int xinc, yinc, zinc;
  int n=0;

  for(zinc=-2;zinc<=2;zinc++)
  {
    zi = patchData[pid].zi + zinc;
    if ((zi < 0) || (zi >= zDim))
      continue;
    for(yinc=-2;yinc<=2;yinc++)
    {
      yi = patchData[pid].yi + yinc;
      if ((yi < 0) || (yi >= yDim))
	continue;
      for(xinc=-2;xinc<=2;xinc++)
      {
	if (!((xinc==2) || (yinc==2) || (zinc==2) ||
	    (xinc==-2) || (yinc==-2) || (zinc==-2)))
	  continue;

	xi = patchData[pid].xi + xinc;
	if ((xi < 0) || (xi >= xDim))
	  continue;

	neighbor_ids[n]=this->pid(xi,yi,zi);
	n++;
      }
    }
  }
  return n;
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
  patchData[pid].myPatch = NULL;
}



/***************************************************************************
 * RCS INFORMATION:
 *
 *	$RCSfile: PatchMap.C,v $
 *	$Author: jim $	$Locker:  $		$State: Exp $
 *	$Revision: 1.10 $	$Date: 1996/12/12 08:57:17 $
 *
 ***************************************************************************
 * REVISION HISTORY:
 *
 * $Log: PatchMap.C,v $
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
