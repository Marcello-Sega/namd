/**
***  Copyright (c) 1995, 1996, 1997, 1998, 1999, 2000 by
***  The Board of Trustees of the University of Illinois.
***  All rights reserved.
**/

#include <stddef.h>
#if !defined(WIN32) || defined(__CYGWIN__)
#include <unistd.h>
#endif
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
  nNodesWithPatches = 0;
  int npes = CkNumPes();
  nPatchesOnNode = new int[npes];
  for ( int i=0; i<npes; ++i ) {
    nPatchesOnNode[i] = 0;
  }
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
    int max_computes = 200;
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
  delete [] nPatchesOnNode;
  
}

#undef PACK
#define PACK(type,data) { memcpy(b, &data,sizeof(type)); b += sizeof(type); }
#define PACKN(type,data,cnt) { memcpy(b, data,(cnt)*sizeof(type)); b += (cnt)*sizeof(type); }

int PatchMap::packSize(void)
{
  int i, size = 0;
  size += 11 * sizeof(int) + 6 * sizeof(BigReal);
  size += CkNumPes() * sizeof(int);
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
  PACK(int,nNodesWithPatches);
  PACKN(int,nPatchesOnNode,CkNumPes());
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
#define UNPACKN(type,data,cnt) { memcpy(data, b, (cnt)*sizeof(type)); b += (cnt)*sizeof(type); }

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
  UNPACK(int,nNodesWithPatches);
  UNPACKN(int,nPatchesOnNode,CkNumPes());
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
void PatchMap::homePatchIDList(PatchIDList &pids) {
  pids.resize(0);
  int i;
  for ( i=0; i<nPatches; ++i ) {
    if ( patchData[i].node == CkMyPe() ) {
      pids.add(i);
    }
  }
}

//----------------------------------------------------------------------
void PatchMap::assignNode(PatchID pid, NodeID node) {
  patchData[pid].node=node;
  if ( nPatchesOnNode[node] == 0 ) nNodesWithPatches += 1;
  nPatchesOnNode[node] += 1;
}

//----------------------------------------------------------------------
void PatchMap::newCid(int pid, int cid)
{
  if (patchData[pid].numCids >= patchData[pid].numCidsAllocated)
  { // allocate more
//    NAMD_die("PatchMap::newCid - not enough compute ID's allocated.");
    ComputeID *old = patchData[pid].cids;
    patchData[pid].numCidsAllocated += 200;
    patchData[pid].cids = new int[patchData[pid].numCidsAllocated];
    int i;
    for (i=0; i<patchData[pid].numCids; i++) 
    	patchData[pid].cids[i] = old[i];
    for (i=patchData[pid].numCids; i<patchData[pid].numCidsAllocated; i++) 
	patchData[pid].cids[i] = -1;
    delete old;
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

