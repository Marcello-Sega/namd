/*****************************************************************************
 * $Source: /home/cvs/namd/cvsroot/namd2/src/TorusLB.C,v $
 * $Author: bhatele $
 * $Date: 2007/11/05 20:23:48 $
 * $Revision: 1.4 $
 *****************************************************************************/
 
/** \file TorusLB.C
 *  Author: Abhinav S Bhatele
 *  Date Created: June 05th, 2007 
 *
 *  Replacement for AlgSeven.C
 */

#include "TorusLB.h"
#include "ProxyMgr.h"
#define EXPAND_INNER_BRICK 2

TorusLB::TorusLB(computeInfo *cs, patchInfo *pas, processorInfo *pes, int ncs, 
int npas, int npes) : RefineTorusLB(cs, pas, pes, ncs, npas, npes, 0)
{
  strategyName = "TorusLB";
  strategy();
  binaryRefine();
  computeAverage();
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
  overLoad = 1.1;

  for(int I=0; I<numComputes; I++) {

  c = (computeInfo *) computePairHeap->deleteMax();
  if ( ! c ) c = (computeInfo *) computeSelfHeap->deleteMax(); 

  if(c->processor != -1) continue; // go to the next compute
  if(!c) CkAbort("TorusLB: Compute Heap empty!\n");

  for(int j=0; j<6; j++) {
    bestPe[j] = 0;
    goodPe[j] = 0;
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
  if((p = bestPe[3])
  || (p = bestPe[4])
  || (p = bestPe[5])
#if USE_TOPOMAP
  || (p = goodPe[3])
  || (p = goodPe[4])
  || (p = goodPe[5])
#endif
  || (p = bestPe[1])
  || (p = bestPe[2])
  || (p = bestPe[0])
#if USE_TOPOMAP
  || (p = goodPe[1])
  || (p = goodPe[2])
  || (p = goodPe[0])
#endif
  ) {
    assign(c, p);
    continue;
  }
 
#if USE_TOPOMAP
  // If no processor found, go through the whole list in a topological fashion
  // first try the inner brick
  int found = 0;
  int p1, p2, pe, x1, x2, xm, xM, y1, y2, ym, yM, z1, z2, zm, zM;
  int dimX, dimY, dimZ;
  double minLoad;
  p1 = patches[c->patch1].processor;
  p2 = patches[c->patch2].processor;

  tmgr.rankToCoordinates(p1, x1, y1, z1);
  tmgr.rankToCoordinates(p2, x2, y2, z2);
  dimX = tmgr.getDimX();
  dimY = tmgr.getDimY();
  dimZ = tmgr.getDimZ();

  brickDim(x1, x2, dimX, xm, xM);
  brickDim(y1, y2, dimY, ym, yM);
  brickDim(z1, z2, dimZ, zm, zM);

  // to expand the inner brick by some hops
#if 0
  if(xm>=EXPAND_INNER_BRICK) xm=xm-EXPAND_INNER_BRICK; else xm=0;
  if(ym>=EXPAND_INNER_BRICK) ym=ym-EXPAND_INNER_BRICK; else ym=0;
  if(zm>=EXPAND_INNER_BRICK) zm=zm-EXPAND_INNER_BRICK; else zm=0;

  if(xM<tmgr.getDimX()-EXPAND_INNER_BRICK) xM=xM+EXPAND_INNER_BRICK; else xM=tmgr.getDimX()-1;
  if(yM<tmgr.getDimY()-EXPAND_INNER_BRICK) yM=yM+EXPAND_INNER_BRICK; else yM=tmgr.getDimY()-1;
  if(zM<tmgr.getDimZ()-EXPAND_INNER_BRICK) zM=zM+EXPAND_INNER_BRICK; else zM=tmgr.getDimZ()-1;
#endif

  // first go over the processors inside the brick and choose the least 
  // overloaded one
  p = 0; minp = 0;
  minLoad = overLoad * averageLoad;
  for(int i=xm; i<=xM; i++)
    for(int j=ym; j<=yM; j++)
      for(int k=zm; k<=zM; k++)
      {
        pe = tmgr.coordinatesToRank(i%dimX, j%dimY, k%dimZ);
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
    for(int i=xM+1; i<xm+dimX; i++)
      for(int j=yM+1; j<ym+dimY; j++)
        for(int k=zM+1; k<zm+dimZ; k++)
        {
          pe = tmgr.coordinatesToRank(i%dimX, j%dimY, k%dimZ);
          p = &processors[pe];
          if(c->load + p->load < minLoad) { 
            minp = p;
	    found = 1; break;
          }
        }
  }
  
  if(found == 0)
    CkAbort("TorusLB: No receiver found\n");
  else  
    assign(c, minp);
#else
  int found = 0;
  if(!minp) {
    heapIterator nextp;
    processorInfo *p = (processorInfo *)(pes->iterator((heapIterator *) &nextp));
    while (p) {
      selectPes(p, c);
      p = (processorInfo *)(pes->next(&nextp));
    }
    p = 0;
    if((p = bestPe[5])
    || (p = bestPe[4])
    || (p = bestPe[3])
    || (p = bestPe[2])
    || (p = bestPe[1])
    || (p = bestPe[0])) {
      assign(c, p);
      found = 1;
      continue;
    }
  }

  if(found == 0)
     CkAbort("TorusLB: No receiver found\n");
#endif
 
  } // end of computes for-loop

  computeAverage();
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

#if USE_TOPOMAP
  int x, y, z;
  int p1, p2, pe, x1, x2, xm, xM, y1, y2, ym, yM, z1, z2, zm, zM;
  int dimX, dimY, dimZ;
  double minLoad;
  p1 = patches[c->patch1].processor;
  p2 = patches[c->patch2].processor;

  tmgr.rankToCoordinates(p1, x1, y1, z1);
  tmgr.rankToCoordinates(p2, x2, y2, z2);
  dimX = tmgr.getDimX();
  dimY = tmgr.getDimY();
  dimZ = tmgr.getDimZ();

  brickDim(x1, x2, dimX, xm, xM);
  brickDim(y1, y2, dimY, ym, yM);
  brickDim(z1, z2, dimZ, zm, zM);
#endif

  if(p->load + c->load < overLoad * averageLoad) {
#if USE_TOPOMAP
    tmgr.rankToCoordinates(p->Id, x, y, z);
    int wB = withinBrick(x, y, z, xm, xM, dimX, ym, yM, dimY, zm, zM, dimZ);
    if(wB)  {
#endif
      processorInfo* &newp = bestPe[index];
      if (!(newp) || p->load < newp->load )
        newp = p;
#if USE_TOPOMAP
    }
    else {
      processorInfo* &newp = goodPe[index];
      if (!(newp) || p->load < newp->load )
        newp = p;
    }
#endif
  }
}


