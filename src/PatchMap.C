/***************************************************************************/
/*                                                                         */
/*              (C) Copyright 1996 The Board of Trustees of the            */
/*                          University of Illinois                         */
/*                           All Rights Reserved                           */
/*									   */
/***************************************************************************/

/***************************************************************************
 * RCS INFORMATION:
 *
 *	$RCSfile: PatchMap.C,v $
 *	$Author: brunner $	$Locker:  $		$State: Exp $
 *	$Revision: 1.1 $	$Date: 1996/08/16 20:43:53 $
 *
 ***************************************************************************
 * DESCRIPTION:
 *
 ***************************************************************************
 * REVISION HISTORY:
 *
 * $Log: PatchMap.C,v $
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

#include <stddef.h>
#include <unistd.h>
#include "chare.h"
#include "ckdefs.h"
#include "c++interface.h"

#include "PatchMap.h"

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
  return ((zIndex*zDim)+yIndex)*yDim + xIndex;
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
  return (pid / xDim) / yDim;
}

//----------------------------------------------------------------------
int PatchMap::node(int pid)
{
  return patchData[pid].node;
}

//----------------------------------------------------------------------
Position PatchMap::minX(int pid)
{
  return patchData[pid].x0;
}

//----------------------------------------------------------------------
Position PatchMap::maxX(int pid)
{
  return patchData[pid].x1;
}

//----------------------------------------------------------------------
Position PatchMap::minY(int pid)
{
  return patchData[pid].y0;
}

//----------------------------------------------------------------------
Position PatchMap::maxY(int pid)
{
  return patchData[pid].y1;
}

//----------------------------------------------------------------------
Position PatchMap::minZ(int pid)
{
  return patchData[pid].z0;
}

//----------------------------------------------------------------------
Position PatchMap::maxZ(int pid)
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
int PatchMap::allocatePids(int ixDim, int iyDim, int izDim)
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
    return -1;

  for(i=0;i<nPatches;i++)
  {
    patchData[i].numCids=0;
    patchData[i].xi = xIndex(i);
    patchData[i].yi = yIndex(i);
    patchData[i].zi = zIndex(i);
    patchData[i].cids=NULL;
    patchData[i].numCidsAllocated=0;
  }

  return 0;
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
			  Position x0, Position y0, Position z0,
			  Position x1, Position y1, Position z1)
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
int PatchMap::newCid(int pid, int cid)
{
  if (patchData[pid].numCids < patchData[pid].numCidsAllocated)
  {
    patchData[pid].cids[patchData[pid].numCids]=cid;
    patchData[pid].numCids++;
  } else return -1;
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
