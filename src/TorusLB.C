/*****************************************************************************
 * $Source: /home/cvs/namd/cvsroot/namd2/src/TorusLB.C,v $
 * $Author: bhatele $
 * $Date: 2008/08/27 02:35:17 $
 * $Revision: 1.13 $
 *****************************************************************************/
 
/** \file TorusLB.C
 *  Author: Abhinav S Bhatele
 *  Date Created: June 05th, 2007 
 *
 *  Replacement for AlgSeven.C
 */

#include "TorusLB.h"
#include "ProxyMgr.h"
#define SHRINK_INNER_BRICK 1

TorusLB::TorusLB(computeInfo *cs, patchInfo *pas, processorInfo *pes, int ncs, 
int npas, int npes) : RefineTorusLB(cs, pas, pes, ncs, npas, npes, 0)
{
  strategyName = "TorusLB";
  strategy();
  binaryRefine();
  printLoads();
  // CREATE THE SPANNING TREE IN THE LOAD BALANCER
  //if(proxySendSpanning || proxyRecvSpanning)
  //  createSpanningTree();
}

TorusLB::~TorusLB() { }

void TorusLB::strategy() {
  // compute the average load by (compute load + background load) / numPesAvailable
  computeAverage();
  // two heaps of self and pair computes
  makeTwoHeaps();

  computeInfo *c;
  processorInfo *p, *minp;
  Iterator nextP;
  overLoad = 1.2;

  for(int I=0; I<numComputes; I++) {

  c = (computeInfo *) computePairHeap->deleteMax();
  if ( ! c ) c = (computeInfo *) computeSelfHeap->deleteMax(); 

  if(c->processor != -1) continue; // go to the next compute
  if(!c) CkAbort("TorusLB: Compute Heap empty!\n");

  for(int j=0; j<6; j++) {
    bestPe[j] = 0;
    goodPe[j] = 0;
    badPe[j] = 0;
  }

  // Look at the processors which have the compute's patches first
  p = &processors[patches[c->patch1].processor];	// patch 1
  selectPes(p, c);
  p = &processors[patches[c->patch2].processor];	// patch 2
  selectPes(p, c); 

  // Try the processors which have the patches' proxies
  p = (processorInfo *)(patches[c->patch1].proxiesOn.iterator((Iterator *)&nextP));
  while(p) {						// patch 1
    selectPes(p, c);
    p = (processorInfo *)(patches[c->patch1].proxiesOn.next((Iterator *)&nextP));
  } 

  p = (processorInfo *)(patches[c->patch2].proxiesOn.iterator((Iterator *)&nextP));
  while(p) {						// patch 2
    selectPes(p, c);
    p = (processorInfo *)(patches[c->patch2].proxiesOn.next((Iterator *)&nextP));
  }

  // see if we have found a processor to place the compute on
  p = 0;
  if((p = bestPe[5])
#if USE_TOPOMAP
  || (p = goodPe[5])
#endif
  || (p = bestPe[4])
#if USE_TOPOMAP
  || (p = goodPe[4])
#endif
  || (p = bestPe[3])
#if USE_TOPOMAP
  || (p = goodPe[3])
#endif
  || (p = bestPe[1])
#if USE_TOPOMAP
  || (p = goodPe[1])
#endif
  || (p = bestPe[2])
#if USE_TOPOMAP
  || (p = goodPe[2])
#endif
  || (p = bestPe[0])
#if USE_TOPOMAP
  || (p = goodPe[0])
#endif
  ) {
    assign(c, p);
    continue;
  }
 
  int found = 0;
#if USE_TOPOMAP
  // If no processor found, go through the whole list in a topological fashion
  // first try the inner brick
  int p1, p2, pe, x1, x2, xm, xM, y1, y2, ym, yM, z1, z2, zm, zM, t1, t2;
  int dimNX, dimNY, dimNZ, dimNT;
  double minLoad;
  p1 = patches[c->patch1].processor;
  p2 = patches[c->patch2].processor;

  tmgr.rankToCoordinates(p1, x1, y1, z1, t1);
  tmgr.rankToCoordinates(p2, x2, y2, z2, t2);
  dimNX = tmgr.getDimNX();
  dimNY = tmgr.getDimNY();
  dimNZ = tmgr.getDimNZ();
  dimNT = tmgr.getDimNT();

  brickDim(x1, x2, dimNX, xm, xM);
  brickDim(y1, y2, dimNY, ym, yM);
  brickDim(z1, z2, dimNZ, zm, zM);

  // to shrink the inner brick by some hops
#if 0
  xm=xm+SHRINK_INNER_BRICK;
  ym=ym+SHRINK_INNER_BRICK;
  zm=zm+SHRINK_INNER_BRICK;

  xM=xM-SHRINK_INNER_BRICK;
  yM=yM-SHRINK_INNER_BRICK;
  zM=zM-SHRINK_INNER_BRICK;
#endif

  // first go over the processors inside the brick and choose the least 
  // overloaded one
  p = 0; minp = 0;
  minLoad = overLoad * averageLoad;
  for(int i=xm; i<=xM; i++)
    for(int j=ym; j<=yM; j++)
      for(int k=zm; k<=zM; k++)
        for(int l=0; l<dimNT; l++)
        {
          pe = tmgr.coordinatesToRank(i%dimNX, j%dimNY, k%dimNZ, l);
          p = &processors[pe];
          if(c->load + p->load < minLoad) { 
            minLoad = c->load + p->load;
            minp = p;
	    found = 1;
          }
        }

  // if no success, then go through the remaining torus (outer brick)
  // and pick the first underloaded one
  minLoad = overLoad * averageLoad;
  if(found == 0) {
    p = 0; minp = 0;
    for(int i=xM+1; i<xm+dimNX; i++)
      for(int j=0; j<dimNY; j++)
	for(int k=0; k<dimNZ; k++)
          for(int l=0; l<dimNT; l++)
          {
            pe = tmgr.coordinatesToRank(i%dimNX, j%dimNY, k%dimNZ, l);
            p = &processors[pe];
            if(c->load + p->load < minLoad) { 
              minp = p;
	      found = 1; break;
            }
          }
  }

  if(found == 0) {
    for(int j=yM+1; j<ym+dimNY; j++)
      for(int i=xm; i<=xM; i++)
	for(int k=0; k<dimNZ; k++)
          for(int l=0; l<dimNT; l++)
          {
            pe = tmgr.coordinatesToRank(i%dimNX, j%dimNY, k%dimNZ, l);
            p = &processors[pe];
            if(c->load + p->load < minLoad) { 
              minp = p;
	      found = 1; break;
            }
          }
  }

  if(found == 0) {
    for(int k=zM+1; k<zm+dimNZ; k++)
      for(int i=xm; i<=xM; i++)
        for(int j=ym; j<=yM; j++)
          for(int l=0; l<dimNT; l++)
          {
            pe = tmgr.coordinatesToRank(i%dimNX, j%dimNY, k%dimNZ, l);
            p = &processors[pe];
            if(c->load + p->load < minLoad) { 
              minp = p;
	      found = 1; break;
            }
          }
  }

  
  if(found == 1) {
    assign(c, minp);
    continue;
  }

#endif /* USE_TOPOMAP */

  if(found == 0) {
    heapIterator nextp;
    processorInfo *p = (processorInfo *)(pes->iterator((heapIterator *) &nextp));
    while (p) {
      selectPes(p, c);
      p = (processorInfo *)(pes->next(&nextp));
    }
    p = 0;
    if((p = bestPe[5])
#if USE_TOPOMAP
    || (p = goodPe[5])
#endif
    || (p = bestPe[4])
#if USE_TOPOMAP
    || (p = goodPe[4])
#endif
    || (p = bestPe[3])
#if USE_TOPOMAP
    || (p = goodPe[3])
#endif
    || (p = bestPe[1])
#if USE_TOPOMAP
    || (p = goodPe[1])
#endif
    || (p = bestPe[2])
#if USE_TOPOMAP
    || (p = goodPe[2])
#endif
    || (p = bestPe[0])
#if USE_TOPOMAP
    || (p = goodPe[0])
#endif
    ) {
      assign(c, p);
      found = 1;
      continue;
    }
  }

  if(found == 0) {
    p = 0;
    if((p = badPe[5])
    || (p = badPe[4])
    || (p = badPe[3])
    || (p = badPe[1])
    || (p = badPe[2])
    || (p = badPe[0])) {
      assign(c, p);
      found = 1;
      continue;
    }
  }

  if(found == 0) {
     CkPrintf("TorusLB: No receiver found average %f overload %f\n", averageLoad, overLoad);
     CkAbort("TorusLB: No receiver found\n");
  }
 
  } // end of computes for-loop

  printLoads();
}

void TorusLB::selectPes(processorInfo *p, computeInfo *c) {
  if(p->available == CmiFalse)
    return;

  // find the position in bestPe/goodPe to place this pair
  // HP HP HP HP HP HP
  // 02 11 20 01 10 00
  //  5  4  3  2  1  0
  int numPatches, numProxies, badForComm, index;
  numAvailable(c, p, &numPatches, &numProxies, &badForComm); 
  index = ((numPatches==2) ? (numPatches+1) : numPatches) + (numProxies * 2 + 1);

  if(numPatches==0 && numProxies==1)
    index--;
  if(numProxies==0)
    index--; 
  /*if(numPatches == 2 || numProxies == 2) index = 5;
  if(numPatches == 1 && numProxies == 1) index = 5;
  else if(numPatches == 1 || numProxies == 1) index = 4;
  else index = 3;*/

#if USE_TOPOMAP
  int x, y, z, t;
  int p1, p2, pe, x1, x2, xm, xM, y1, y2, ym, yM, z1, z2, zm, zM, t1, t2;
  int dimNX, dimNY, dimNZ, dimNT;
  double minLoad;
  p1 = patches[c->patch1].processor;
  p2 = patches[c->patch2].processor;

  tmgr.rankToCoordinates(p1, x1, y1, z1, t1);
  tmgr.rankToCoordinates(p2, x2, y2, z2, t2);
  dimNX = tmgr.getDimNX();
  dimNY = tmgr.getDimNY();
  dimNZ = tmgr.getDimNZ();
  dimNT = tmgr.getDimNT();

  brickDim(x1, x2, dimNX, xm, xM);
  brickDim(y1, y2, dimNY, ym, yM);
  brickDim(z1, z2, dimNZ, zm, zM);
#endif

  if(p->load + c->load < overLoad * averageLoad) {
#if USE_TOPOMAP
    tmgr.rankToCoordinates(p->Id, x, y, z, t);
    int wB = withinBrick(x, y, z, xm, xM, dimNX, ym, yM, dimNY, zm, zM, dimNZ);
    if(wB) {
#endif
      processorInfo* &newp = bestPe[index];
      if (!(newp) || p->load < newp->load )
        newp = p;
#if USE_TOPOMAP
    }
    else {
      processorInfo* &newp = goodPe[index];
      if (!(newp) /*|| p->load < newp->load*/ )
        newp = p;
      else {
        if(tmgr.getHopsBetweenRanks(newp->Id, p1) + tmgr.getHopsBetweenRanks(newp->Id, p2) > tmgr.getHopsBetweenRanks(p->Id, p1) + tmgr.getHopsBetweenRanks(p->Id, p2))
          newp = p;
      }
    }
#endif
  }
  else {
    processorInfo* &newp = badPe[index];
    if (!(newp) || p->load < newp->load )
      newp = p;
  }
}


