/**
***  Copyright (c) 1995, 1996, 1997, 1998, 1999, 2000 by
***  The Board of Trustees of the University of Illinois.
***  All rights reserved.
**/

#include <iostream.h>
#include "common.h"
#include "InfoStream.h"
#include "AlgRecBisection.h"

//
//   load balancing computes using recursion of bisection algorithm
//
//
AlgRecBisection::AlgRecBisection(computeInfo *computeArray, patchInfo *patchArray, 
	   processorInfo *processorArray, int nComps, 
	   int nPatches, int nPes) :
  Rebalancer(computeArray, patchArray, 
	     processorArray, nComps, 
	     nPatches, nPes)
{
strategyName = "Rob";
strategy();
}


void AlgRecBisection::rec_divide(int n, Partition &p)
{
  int i,j,k;
  int pos[3];
  int n1, n2;
  double load1, currentload;
  int mindir, count;
  Partition p1, p2;

  CmiPrintf("rec_divide: partition n:%d count:%d load:%f (%d %d %d, %d %d %d)\n", n, p.count, p.load, p.origin[0], p.origin[1], p.origin[2], p.corner[0], p.corner[1], p.corner[2]);

  if (n==1) {
    partitions[currentp++] = p;
    return;
  }
/*
  if (p.origin.x==p.corner.x && p.origin.y==p.corner.y && p.origin.z==p.corner.z) 
     NAMD_die("AlgRecBisection failed in recursion.\n"); 
*/

  n2 = n/2;
  n1 = n-n2;

  load1 = (1.0*n1/n) * p.load;

  p1 = p;
  p1.refno = ++refno;
  p2 = p;
  p2.refno = ++refno;

  // determine the best division direction
  int maxSpan=-1;
  for (i=XDIR; i<=ZDIR; i++) {
    int myspan = p.corner[i] - p.origin[i];
    if (myspan > maxSpan) {
      mindir = i;
      maxSpan = myspan;
    }
  }

  // other two dimensions
  int dir2 = (mindir+1)%3;
  int dir3 = (mindir+2)%3;

  currentload = 0.0;
  count = 0;
  pos[mindir] = p.origin[mindir];
  for (i=0; i<numComputes; i++) {
    // not belong to this partition
    if (computeLoad[vArray[mindir][i].id].refno != p.refno) continue;
    if (vArray[mindir][i].v<p.origin[mindir]) continue;
    if (vArray[mindir][i].v>p.corner[mindir]) break;

    int cid = vArray[mindir][i].id;	// this compute ID
    // check if this compute is within the partition
    if ( computeLoad[cid].v[dir2] >= p.origin[dir2] &&
	 computeLoad[cid].v[dir2] <= p.corner[dir2] &&
	 computeLoad[cid].v[dir3] >= p.origin[dir3] &&
	 computeLoad[cid].v[dir3] <= p.corner[dir3]  ) {
      // this compute is set to the first partition
      if (currentload <= load1) {
	computeLoad[cid].refno = p1.refno;
        currentload += computeLoad[cid].load;
        count ++;
	pos[mindir] = computeLoad[cid].v[mindir];
      }
      else {	// or the next partition
	computeLoad[cid].refno = p2.refno;
      }
    }
  }
//  CmiPrintf("X:cur:%d, prev:%d load:%f %f\n", cur, prev, currentload, prevload);
  CmiPrintf("DIR:%d %d load:%f\n", mindir, pos[mindir], currentload);


  p1.corner[mindir] = pos[mindir];
  p2.origin[mindir] = pos[mindir];

  p1.load = currentload;
  p1.count = count;
  p2.load = p.load - p1.load;
  p2.count = p.count - p1.count;
  CmiPrintf("p1: n:%d count:%d load:%f\n", n1, p1.count, p1.load);
  CmiPrintf("p2: n:%d count:%d load:%f\n", n2, p2.count, p2.load);
  rec_divide(n1, p1);
  rec_divide(n2, p2);
}

int comp(const void *a, const void *b)
{
  AlgRecBisection::VecArray *va = (AlgRecBisection::VecArray *)a;
  AlgRecBisection::VecArray *vb = (AlgRecBisection::VecArray *)b;
  return va->v - vb->v;
}

void AlgRecBisection::strategy()
{
  int i,j;

  PatchMap *patchMap = PatchMap::Object();

  // create computeLoad and calculate tentative computes coordinates
  computeLoad = new ComputeLoad[numComputes];
  vArray[XDIR] = new VecArray[numComputes];
  vArray[YDIR] = new VecArray[numComputes];
  vArray[ZDIR] = new VecArray[numComputes];

  CmiPrintf("AlgRecBisection: numComputes:%d\n", numComputes);

  for (i=0; i<numComputes; i++) {
    int pid1 = computes[i].patch1;
    int pid2 = computes[i].patch2;
    computeLoad[i].id = i;
    computeLoad[i].v[0] = patchMap->index_a(pid1) + patchMap->index_a(pid2);
    computeLoad[i].v[1] = patchMap->index_b(pid1) + patchMap->index_b(pid2);
    computeLoad[i].v[2] = patchMap->index_c(pid1) + patchMap->index_c(pid2);
//    CmiPrintf("(%d %d %d)", computeLoad[i].x, computeLoad[i].y, computeLoad[i].z);
    computeLoad[i].load = computes[i].load;
    computeLoad[i].refno = 0;

    vArray[XDIR][i].id = vArray[YDIR][i].id = vArray[ZDIR][i].id = i;
    vArray[XDIR][i].v = computeLoad[i].v[0];
    vArray[YDIR][i].v = computeLoad[i].v[1];
    vArray[ZDIR][i].v = computeLoad[i].v[2];
  }
//  CmiPrintf("\n");

  double t = CmiWallTimer();
  qsort(vArray[XDIR], numComputes, sizeof(VecArray), comp);
  qsort(vArray[YDIR], numComputes, sizeof(VecArray), comp);
  qsort(vArray[ZDIR], numComputes, sizeof(VecArray), comp);
  CmiPrintf("qsort time: %f\n", CmiWallTimer() - t);

  npartition = P;
  partitions = new Partition[npartition];

  top_partition.origin[XDIR] = 0;
  top_partition.origin[YDIR] = 0;
  top_partition.origin[ZDIR] = 0;
  top_partition.corner[XDIR] = 2*(patchMap->gridsize_a()-1);
  top_partition.corner[YDIR] = 2*(patchMap->gridsize_b()-1);
  top_partition.corner[ZDIR] = 2*(patchMap->gridsize_c()-1);

  top_partition.refno = 0;
  top_partition.load = 0.0;
  top_partition.count = numComputes;
  for (i=0; i<numComputes; i++) top_partition.load += computes[i].load;

  currentp = 0;
  refno = 0;

  // recursively divide
  rec_divide(npartition, top_partition);

  CmiPrintf("After partitioning: \n");
  for (i=0; i<P; i++) {
    CmiPrintf("[%d] (%d,%d,%d) (%d,%d,%d) load:%f count:%d\n", i, partitions[i].origin[0], partitions[i].origin[1], partitions[i].origin[2], partitions[i].corner[0], partitions[i].corner[1], partitions[i].corner[2], partitions[i].load, partitions[i].count);
  }

  for (i=0; i<numComputes; i++) computes[i].processor=-1;

  // this is for debugging
  int *num = new int[P];
  for (i=0; i<P; i++) num[i] = 0;

  for (i=0; i<numComputes; i++)
  {
    for (j=0; j<P; j++)
      if (computeLoad[i].refno == partitions[j].refno) 
        { computes[computeLoad[i].id].processor = j; num[j]++; break; }
  }

  for (i=0; i<P; i++)
    if (num[i] != partitions[i].count) 
      NAMD_die("AlgRecBisection: Compute counts don't agree!\n");

  for (i=0; i<numComputes; i++) {
    if (computes[i].processor == -1) NAMD_bug("AlgRecBisection failure!\n");
//    CmiPrintf("[%d] old:%d new:%d\n", i, computes[i].oldProcessor, computes[i].processor);
  }


  delete [] computeLoad;
  for (i=0; i<3; i++) delete [] vArray[i];
  delete [] partitions;

  CmiPrintf("AlgRecBisection finished\n");
  CmiPrintf("AlgRecBisection finish time: %f\n", CmiWallTimer() - t);

  // printLoads();
}


