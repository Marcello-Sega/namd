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

#include "charm++.h"

#include "PatchMgr.h"
#include "PatchMap.inl"
#include "Patch.h"
#include "Lattice.h"
#include "HomePatchList.h"
#include "AtomMap.h"

//#define DEBUGM
#define MIN_DEBUG_LEVEL 5
#include "Debug.h"

// Safe singleton creation
PatchMap *PatchMap::Instance() {
  if (CpvAccess(PatchMap_instance) == 0) {
     CpvAccess(PatchMap_instance) = new PatchMap;
  }
  return(CpvAccess(PatchMap_instance));
}


PatchMap::PatchMap(void)
{
  nPatches = 0;
  patchData = NULL;
  aDim = bDim = cDim = 0;
  aPeriodic = bPeriodic = cPeriodic = 0;
  aMaxIndex = bMaxIndex = cMaxIndex = 0;
}

void PatchMap::initialize(ScaledPosition xmin, ScaledPosition xmax,
                              Lattice lattice, BigReal patchSize)
{
  aPeriodic = lattice.a_p();
  bPeriodic = lattice.b_p();
  cPeriodic = lattice.c_p();

  if ( aPeriodic ) {
    BigReal sysDim = lattice.a_r().unit() * lattice.a();
    aDim = (int)(sysDim / patchSize);
  } else {
    BigReal sysDim = xmax.x - xmin.x;
    aDim = (int)(sysDim / patchSize);
    if ((aDim * patchSize) < sysDim) aDim++;
  }

  if ( bPeriodic ) {
    BigReal sysDim = lattice.b_r().unit() * lattice.b();
    bDim = (int)(sysDim / patchSize);
  } else {
    BigReal sysDim = xmax.y - xmin.y;
    bDim = (int)(sysDim / patchSize);
    if ((bDim * patchSize) < sysDim) bDim++;
  }

  if ( cPeriodic ) {
    BigReal sysDim = lattice.c_r().unit() * lattice.c();
    cDim = (int)(sysDim / patchSize);
  } else {
    BigReal sysDim = xmax.z - xmin.z;
    cDim = (int)(sysDim / patchSize);
    if ((cDim * patchSize) < sysDim) cDim++;
  }

  if ( aDim < 0 || bDim < 0 || cDim < 0 ) {
    NAMD_die("Bug in PatchMap::initialize - negative grid dimension.");
  }

  if ( aDim == 0 ) aDim = 1;
  if ( bDim == 0 ) bDim = 1;
  if ( cDim == 0 ) cDim = 1;

  if ( aPeriodic && aDim < 2 ) aDim = 2;
  if ( bPeriodic && bDim < 2 ) bDim = 2;
  if ( cPeriodic && cDim < 2 ) cDim = 2;

  iout << iINFO << "PATCH GRID IS ";
  iout << aDim;
  if ( aPeriodic ) iout << " (PERIODIC)";
  iout << " BY ";
  iout << bDim;
  if ( bPeriodic ) iout << " (PERIODIC)";
  iout << " BY ";
  iout << cDim;
  if ( cPeriodic ) iout << " (PERIODIC)";
  iout << "\n" << endi;

  aMaxIndex = ( ! aPeriodic || aDim == 2 ) ? 10000 : aDim;
  bMaxIndex = ( ! bPeriodic || bDim == 2 ) ? 10000 : bDim;
  cMaxIndex = ( ! cPeriodic || cDim == 2 ) ? 10000 : cDim;

  aLength = aPeriodic ? 1.0 : aDim * patchSize;
  bLength = bPeriodic ? 1.0 : bDim * patchSize;
  cLength = cPeriodic ? 1.0 : cDim * patchSize;

  aOrigin = aPeriodic ? -0.5 : 0.5 * (xmin.x + xmax.x - aLength);
  bOrigin = bPeriodic ? -0.5 : 0.5 * (xmin.y + xmax.y - bLength);
  cOrigin = cPeriodic ? -0.5 : 0.5 * (xmin.z + xmax.z - cLength);

  nPatches=aDim*bDim*cDim;
  patchData = new PatchData[nPatches];

  for(int i=0; i<nPatches; ++i)
  {
    PatchData &p = patchData[i];
    p.numCids = 0;
    p.aIndex = index_a(i);
    p.bIndex = index_b(i);
    p.cIndex = index_c(i);
    p.myPatch = 0;
    p.myHomePatch = 0;
    p.aMin = ((float)p.aIndex/(float)aDim) * aLength + aOrigin;
    p.bMin = ((float)p.bIndex/(float)bDim) * bLength + bOrigin;
    p.cMin = ((float)p.cIndex/(float)cDim) * cLength + cOrigin;
    p.aMax = ((float)(p.aIndex+1)/(float)aDim) * aLength + aOrigin;
    p.bMax = ((float)(p.bIndex+1)/(float)bDim) * bLength + bOrigin;
    p.cMax = ((float)(p.cIndex+1)/(float)cDim) * cLength + cOrigin;
    p.numCids = 0;
    int max_computes = 100;
    p.cids = new int[max_computes];
    for ( int j = 0; j < max_computes; ++j ) p.cids[j] = -1;
    p.numCidsAllocated = max_computes;
  }

}

void PatchMap::checkMap(void)
{
  int patchCount=0;
  for (int i=0; i<nPatches; i++) {
    if (patchData[i].myPatch) {
      patchCount++;
      if ( patchData[i].myPatch->getPatchID() != i) {
	DebugM(4, "patchID("<<patchData[i].myPatch->getPatchID()
	  <<") != patchID(" 
	  <<i<<")\n");
      }
    }
  }
  DebugM(4, "Patch Count = " <<patchCount<<"\n");
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
#define PACK(type,data) { memcpy(b, &data,sizeof(type)); b += sizeof(type); }

int PatchMap::packSize(void)
{
  int i, size = 0;
  size += 11 * sizeof(int) + 6 * sizeof(BigReal);
  for(i=0;i<nPatches;++i)
  {
    size += sizeof(PatchData);
    size += patchData[i].numCidsAllocated * sizeof(ComputeID);
  }
  return size;
}

void PatchMap::pack (char *buffer)
{
  DebugM(4,"Packing PatchMap on node " << CkMyPe() << endl);
  int i,j;

  // fill in the data
  char *b = buffer;
  PACK(int,nPatches);
  DebugM(3,"nPatches = " << nPatches << endl);
  PACK(int,aDim); PACK(int,bDim); PACK(int,cDim);
  PACK(int,aPeriodic); PACK(int,bPeriodic); PACK(int,cPeriodic);
  PACK(int,aMaxIndex); PACK(int,bMaxIndex); PACK(int,cMaxIndex);
  PACK(BigReal,aOrigin); PACK(BigReal,bOrigin); PACK(BigReal,cOrigin);
  PACK(BigReal,aLength); PACK(BigReal,bLength); PACK(BigReal,cLength);
  for(i=0;i<nPatches;++i)
  {
    DebugM(3,"Packing Patch " << i << " is on node " << patchData[i].node << 
	" with " << patchData[i].numCidsAllocated << " allocated.\n");
    PACK(PatchData,patchData[i]);
    for(j=0;j<patchData[i].numCidsAllocated;++j)
      PACK(ComputeID,patchData[i].cids[j]);
  }
  DebugM(3,buffer + size - b << " == 0 ?" << endl);
}

#undef UNPACK
#define UNPACK(type,data) { memcpy(&data, b, sizeof(type)); b += sizeof(type); }

void PatchMap::unpack (char *ptr)
{
  DebugM(4,"Unpacking PatchMap on node " << CkMyPe() << endl);
  int i,j;
  char *b = (char*)ptr;
  UNPACK(int,nPatches);
  DebugM(3,"nPatches = " << nPatches << endl);
  UNPACK(int,aDim); UNPACK(int,bDim); UNPACK(int,cDim);
  UNPACK(int,aPeriodic); UNPACK(int,bPeriodic); UNPACK(int,cPeriodic);
  UNPACK(int,aMaxIndex); UNPACK(int,bMaxIndex); UNPACK(int,cMaxIndex);
  UNPACK(BigReal,aOrigin); UNPACK(BigReal,bOrigin); UNPACK(BigReal,cOrigin);
  UNPACK(BigReal,aLength); UNPACK(BigReal,bLength); UNPACK(BigReal,cLength);
  patchData = new PatchData[nPatches];
  for(i=0;i<nPatches;++i)
  {
    UNPACK(PatchData,patchData[i]);
    if (CkMyPe()) {
      patchData[i].myPatch = 0;
      patchData[i].myHomePatch = 0;
    }
    DebugM(3,"Unpacking Patch " << i << " is on node " << patchData[i].node << 
	" with " << patchData[i].numCidsAllocated << " allocated.\n");
    patchData[i].cids = new ComputeID[patchData[i].numCidsAllocated];
    for(j=0;j<patchData[i].numCidsAllocated;++j)
      UNPACK(ComputeID,patchData[i].cids[j]);
  }
}

//----------------------------------------------------------------------
int PatchMap::numHomePatches(void)
{
  return CpvAccess(PatchMap_patchMgr)->homePatches.size();
}

//----------------------------------------------------------------------
HomePatchList *PatchMap::homePatchList() {
  return &(CpvAccess(PatchMap_patchMgr)->homePatches);
}

//----------------------------------------------------------------------
void PatchMap::assignNode(PatchID pid, NodeID node) {
  patchData[pid].node=node;
}

//----------------------------------------------------------------------
void PatchMap::newCid(int pid, int cid)
{
  if (patchData[pid].numCids >= patchData[pid].numCidsAllocated)
  { // allocate more
    NAMD_die("PatchMap::newCid - not enough compute ID's allocated.");
  }
  patchData[pid].cids[patchData[pid].numCids]=cid;
  patchData[pid].numCids++;
}

//----------------------------------------------------------------------
int PatchMap::oneAwayNeighbors(int pid, PatchID *neighbor_ids, int *transform_ids)
{
  int xi, yi, zi;
  int xinc, yinc, zinc;
  int n=0;

  for(zinc=-1;zinc<=1;zinc++)
  {
    zi = patchData[pid].cIndex + zinc;
    if ((zi < 0) || (zi >= cDim))
      if ( ! cPeriodic ) continue;
    for(yinc=-1;yinc<=1;yinc++)
    {
      yi = patchData[pid].bIndex + yinc;
      if ((yi < 0) || (yi >= bDim))
	if ( ! bPeriodic ) continue;
      for(xinc=-1;xinc<=1;xinc++)
      {
	if ((xinc==0) && (yinc==0) && (zinc==0))
	  continue;

	xi = patchData[pid].aIndex + xinc;
	if ((xi < 0) || (xi >= aDim))
	  if ( ! aPeriodic ) continue;

	if (neighbor_ids)
	  neighbor_ids[n]=this->pid(xi,yi,zi);
	if ( transform_ids )
	{
	  int xt = 0; if ( xi < 0 ) xt = -1; if ( xi >= aDim ) xt = 1;
	  int yt = 0; if ( yi < 0 ) yt = -1; if ( yi >= bDim ) yt = 1;
	  int zt = 0; if ( zi < 0 ) zt = -1; if ( zi >= cDim ) zt = 1;
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
    zi = patchData[pid].cIndex + zinc;
    if ((zi < 0) || (zi >= cDim))
      if ( ! cPeriodic ) continue;
    for(yinc=-2;yinc<=2;yinc++)
    {
      yi = patchData[pid].bIndex + yinc;
      if ((yi < 0) || (yi >= bDim))
	if ( ! bPeriodic ) continue;
      for(xinc=-2;xinc<=2;xinc++)
      {
	if (!((xinc==2) || (yinc==2) || (zinc==2) ||
	    (xinc==-2) || (yinc==-2) || (zinc==-2)))
	  continue;

	xi = patchData[pid].aIndex + xinc;
	if ((xi < 0) || (xi >= aDim))
	  if ( ! aPeriodic ) continue;

	neighbor_ids[n]=this->pid(xi,yi,zi);
	if ( transform_ids )
	{
	  int xt = 0; if ( xi < 0 ) xt = 1; if ( xi >= aDim ) xt = -1;
	  int yt = 0; if ( yi < 0 ) yt = 1; if ( yi >= bDim ) yt = -1;
	  int zt = 0; if ( zi < 0 ) zt = 1; if ( zi >= cDim ) zt = -1;
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
int PatchMap::upstreamNeighbors(int pid, PatchID *neighbor_ids, int *transform_ids)
{
  int xi, yi, zi;
  int xinc, yinc, zinc;
  int n=0;

  for(zinc=0;zinc<=1;zinc++)
  {
    zi = patchData[pid].cIndex + zinc;
    if ((zi < 0) || (zi >= cDim))
      if ( ! cPeriodic ) continue;
    for(yinc=0;yinc<=1;yinc++)
    {
      yi = patchData[pid].bIndex + yinc;
      if ((yi < 0) || (yi >= bDim))
	if ( ! bPeriodic ) continue;
      for(xinc=0;xinc<=1;xinc++)
      {
	if ((xinc==0) && (yinc==0) && (zinc==0))
	  continue;

	xi = patchData[pid].aIndex + xinc;
	if ((xi < 0) || (xi >= aDim))
	  if ( ! aPeriodic ) continue;

	if (neighbor_ids)
	  neighbor_ids[n]=this->pid(xi,yi,zi);
	if ( transform_ids )
	{
	  int xt = 0; if ( xi < 0 ) xt = -1; if ( xi >= aDim ) xt = 1;
	  int yt = 0; if ( yi < 0 ) yt = -1; if ( yi >= bDim ) yt = 1;
	  int zt = 0; if ( zi < 0 ) zt = -1; if ( zi >= cDim ) zt = 1;
	  transform_ids[n] = Lattice::index(xt,yt,zt);
	}
	n++;
      }
    }
  }
  DebugM(3,"Patch " << pid << " has " << n << " upstream neighbors.\n");
  return n;
}

//----------------------------------------------------------------------
int PatchMap::downstreamNeighbors(int pid, PatchID *neighbor_ids,
				  int *transform_ids)
{
  int xi, yi, zi;
  int xinc, yinc, zinc;
  int n=0;

  for(zinc=-1;zinc<=0;zinc++)
  {
    zi = patchData[pid].cIndex + zinc;
    if ((zi < 0) || (zi >= cDim))
      if ( ! cPeriodic ) continue;
    for(yinc=-1;yinc<=0;yinc++)
    {
      yi = patchData[pid].bIndex + yinc;
      if ((yi < 0) || (yi >= bDim))
	if ( ! bPeriodic ) continue;
      for(xinc=-1;xinc<=0;xinc++)
      {
	if ((xinc==0) && (yinc==0) && (zinc==0))
	  continue;

	xi = patchData[pid].aIndex + xinc;
	if ((xi < 0) || (xi >= aDim))
	  if ( ! aPeriodic ) continue;

	if (neighbor_ids)
	  neighbor_ids[n]=this->pid(xi,yi,zi);
	if ( transform_ids )
	{
	  int xt = 0; if ( xi < 0 ) xt = -1; if ( xi >= aDim ) xt = 1;
	  int yt = 0; if ( yi < 0 ) yt = -1; if ( yi >= bDim ) yt = 1;
	  int zt = 0; if ( zi < 0 ) zt = -1; if ( zi >= cDim ) zt = 1;
	  transform_ids[n] = Lattice::index(xt,yt,zt);
	}
	n++;
      }
    }
  }
  DebugM(3,"Patch " << pid << " has " << n << " upstream neighbors.\n");
  return n;
}

//----------------------------------------------------------------------
void PatchMap::printPatchMap(void)
{
  CkPrintf("---------------------------------------");
  CkPrintf("---------------------------------------\n");

  CkPrintf("nPatches = %d\n",nPatches);
  for(int i=0;i<nPatches;i++)
  {
    CkPrintf("Patch %d:\n",i);
    CkPrintf("  node = %d\n",patchData[i].node);
    CkPrintf("  xi,yi,zi = %d, %d, %d\n",
	    patchData[i].aIndex,patchData[i].bIndex,patchData[i].cIndex);
    CkPrintf("  x0,y0,z0 = %f, %f, %f\n",
	    patchData[i].aMin,patchData[i].bMin,patchData[i].cMin);
    CkPrintf("  x1,y1,z1 = %f, %f, %f\n",
	    patchData[i].aMax,patchData[i].bMax,patchData[i].cMax);
    CkPrintf("  numCids = %d\n",patchData[i].numCids);
    CkPrintf("  numCidsAllocated = %d\n",patchData[i].numCidsAllocated);
    for(int j=0; j < patchData[i].numCids; j++)
    {
      CkPrintf(" %10d ",patchData[i].cids[j]);
      if (!((j+1) % 6))
	CkPrintf("\n");
    }
    CkPrintf("\n---------------------------------------");
    CkPrintf("---------------------------------------\n");
  }

}

//----------------------------------------------------------------------
void PatchMap::registerPatch(PatchID pid, HomePatch *pptr) {
  registerPatch(pid,(Patch*)pptr);
  if (patchData[pid].myHomePatch != 0) {
    iout << iPE << iERRORF 
      << "homePatchID("<<pid<<") is being re-registered!\n" << endi;
  }
  patchData[pid].myHomePatch = pptr;
}

//----------------------------------------------------------------------
void PatchMap::unregisterPatch(PatchID pid, HomePatch *pptr) {
  unregisterPatch(pid,(Patch*)pptr);
  if (pptr == patchData[pid].myHomePatch) {
      DebugM(4, "UnregisterHomePatch("<<pid<<") at " << pptr << "\n");
      patchData[pid].myHomePatch = NULL;
  }
}

//----------------------------------------------------------------------
void PatchMap::registerPatch(PatchID pid, Patch *pptr)
{
  if (patchData[pid].myPatch != 0) {
    iout << iPE << iERRORF 
      << "patchID("<<pid<<") is being re-registered!\n" << endi;
  }
  patchData[pid].myPatch = pptr;
}

//----------------------------------------------------------------------
void PatchMap::unregisterPatch(PatchID pid, Patch *pptr)
{
  if (pptr == patchData[pid].myPatch) {
      DebugM(4, "UnregisterPatch("<<pid<<") at " << pptr << "\n");
      patchData[pid].myPatch = NULL;
  }
}

//----------------------------------------------------------------------
HomePatch *PatchMap::homePatch(PatchID pid)
{
  return patchData[pid].myHomePatch;
}


/***************************************************************************
 * RCS INFORMATION:
 *
 *	$RCSfile: PatchMap.C,v $
 *	$Author: jim $	$Locker:  $		$State: Exp $
 *	$Revision: 1.1022 $	$Date: 1999/09/03 20:46:18 $
 *
 ***************************************************************************
 * REVISION HISTORY:
 *
 * $Log: PatchMap.C,v $
 * Revision 1.1022  1999/09/03 20:46:18  jim
 * Support for non-orthogonal periodic boundary conditions.
 *
 * Revision 1.1021  1999/08/11 16:52:21  jim
 * Make homePatch() method return NULL for proxies rather than casting.
 *
 * Revision 1.1020  1999/05/11 23:56:42  brunner
 * Changes for new charm version
 *
 * Revision 1.1019  1998/07/16 18:52:12  jim
 * Localized common downstream patch optimization.
 *
 * Revision 1.1018  1998/06/24 23:40:33  brunner
 * Added downstreamNeighbors() and LdbCoordinator fixes.  I don't know
 * why Patch.C is different.
 *
 * Revision 1.1017  1998/06/18 14:44:35  jim
 * Eliminated warnings and errors from aCC.
 *
 * Revision 1.1016  1998/03/03 23:05:21  brunner
 * Changed include files for new simplified Charm++ include file structure.
 *
 * Revision 1.1015  1997/11/07 20:17:43  milind
 * Made NAMD to run on shared memory machines.
 *
 * Revision 1.1014  1997/10/06 00:12:33  jim
 * Added PatchMap.inl, sped up cycle-boundary tuple code.
 *
 * Revision 1.1013  1997/09/30 16:57:45  jim
 * Fixed bug dealing with atoms on unknown patches.
 *
 * Revision 1.1012  1997/09/28 22:36:52  jim
 * Modified tuple-based computations to not duplicate calculations and
 * only require "upstream" proxies.
 *
 * Revision 1.1011  1997/04/10 09:14:01  ari
 * Final debugging for compute migration / proxy creation for load balancing.
 * Lots of debug code added, mostly turned off now.
 * Fixed bug in PositionBox when Patch had no dependencies.
 * Eliminated use of cout and misuse of iout in numerous places.
 *                                            Ari & Jim
 *
 * Revision 1.1010  1997/04/06 22:45:08  ari
 * Add priorities to messages.  Mods to help proxies without computes.
 * Added quick enhancement to end of list insertion of ResizeArray(s)
 *
 * Revision 1.1009  1997/04/04 23:34:24  milind
 * Got NAMD2 to run on Origin2000.
 * Included definitions of class static variables in C files.
 * Fixed alignment bugs by using memcpy instead of assignment in
 * pack and unpack.
 *
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
