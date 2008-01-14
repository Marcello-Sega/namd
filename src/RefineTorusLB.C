/*****************************************************************************
 * $Source: /home/cvs/namd/cvsroot/namd2/src/RefineTorusLB.C,v $
 * $Author: bhatele $
 * $Date: 2008/01/14 19:30:41 $
 * $Revision: 1.10 $
 *****************************************************************************/

/** \file RefineTorusLB.C
 *  Author: Abhinav S Bhatele
 *  Date Created: June 12th, 2007 
 *
 *  Replacement for RefineOnly.C
 */

#include "RefineTorusLB.h"
#define EXPAND_INNER_BRICK 2

RefineTorusLB::RefineTorusLB(computeInfo *cs, patchInfo *pas, processorInfo *pes, int ncs, 
int npas, int npes, int flag) : Rebalancer(cs, pas, pes, ncs, npas, npes)
{
  if(flag==1) {
    strategyName = "RefineTorusLB";
    strategy();
    // CREATE THE SPANNING TREE IN THE LOAD BALANCER
#if 0
    if(proxySendSpanning || proxyRecvSpanning) {
    for(int i=0; i<4; i++) {
      decrSTLoad();
      computeAverage();
      createSpanningTree();
      incrSTLoad();
      for(int i=0; i<P; i++)
<     	delete [] processors[i].proxyUsage;
      InitProxyUsage();
      binaryRefine();
      computeAverage();
      printLoads();
      createSpanningTree();
    }
    }
#endif
  }
}

RefineTorusLB::~RefineTorusLB() { }

void RefineTorusLB::strategy() {
  firstAssignInRefine = 0;
  for(int i=0; i<numComputes; i++)
    assign((computeInfo *) &(computes[i]), (processorInfo *) &(processors[computes[i].oldProcessor]));
  firstAssignInRefine = 1;

  binaryRefine();
  
  computeAverage();
  printLoads();
}

void RefineTorusLB::binaryRefine() {
  // CkPrintf("Inside binary refine\n");
  // compute the max and average load
  double avg = computeAverage();
  double max = computeMax();

  double step = 0.01, start = 1.1;
  double dCurLoad = max/avg;
  int curLoad;
  int minLoad = 0;
  int maxLoad = (int)((dCurLoad - start)/step + 1);
  double dMinLoad = minLoad * step + start;
  double dMaxLoad = maxLoad * step + start;

  // check the two limits of the search: 1.1 and dMaxLoad
  int done=0;
  overLoad = dMinLoad;
  if(newRefine())
    done = 1;
  else {
    overLoad = dMaxLoad;
    if(!newRefine()) {
      CkPrintf("Error: Could not refine at max overload\n");
      done = 1;
    }
  } 

  // do a binary search between 1.1 and dMaxLoad until we succeed
  while(!done) {
    //CkPrintf("in multi refine %d %d\n", minLoad, maxLoad);
    if(maxLoad - minLoad <= 1)
      done = 1;
    else {
      curLoad = (maxLoad + minLoad)/2;
      overLoad = curLoad * step + start;
      if(newRefine())
        maxLoad = curLoad;
      else
        minLoad = curLoad;
    }
  }
}

int RefineTorusLB::newRefine() {
  //CkPrintf("Inside new refine %f\n", overLoad);
  int done = 1;
  maxHeap *heavyPes = new maxHeap(P);
  Set *lightPes = new Set();
  processorInfo *donor, *p, *bestP;
  computeInfo *c;
  Iterator nextC, nextP;
  pcpair good;

  // create a heap and set of heavy and light pes respectively
  for(int i=0; i<P; i++) {
    if (processors[i].load > overLoad*averageLoad)
      heavyPes->insert((InfoRecord *) &(processors[i]));
    else
      lightPes->insert((InfoRecord *) &(processors[i]));
  }

  pcpair pcpairarray[12];
     
  for(int j=0; j<6; j++) {
    bestPe[j] = &pcpairarray[j];    // new pcpair();
    goodPe[j] = &pcpairarray[j+6];  // new pcpair();
    //CkPrintf("AllocatE\n");
  }

  while(1) {
    while(donor = (processorInfo*)heavyPes->deleteMax())
      if(donor->computeSet.numElements())
	break;
    
    if(!donor) break;

    for(int j=0; j<6; j++) {
      bestPe[j]->reset();
      goodPe[j]->reset();
      //CkPrintf("AllocatE\n");
    }

    nextC.id = 0;
    c = (computeInfo *)donor->computeSet.iterator((Iterator *)&nextC);

    while(c) {
      // Look at the processors which have the compute's patches first
      p = &processors[patches[c->patch1].processor];        // patch 1
      selectPes(p, c);
      p = &processors[patches[c->patch2].processor];        // patch 2
      selectPes(p, c);

      // Try the processors which have the patches' proxies
      p = (processorInfo *)(patches[c->patch1].proxiesOn.iterator((Iterator *)&nextP));
      while(p) {                                            // patch 1
	selectPes(p, c);
	p = (processorInfo *)(patches[c->patch1].proxiesOn.next((Iterator *)&nextP));
      }
  
      p = (processorInfo *)(patches[c->patch2].proxiesOn.iterator((Iterator *)&nextP));
      while(p) {                                            //patch 2
        selectPes(p, c);
        p = (processorInfo *)(patches[c->patch2].proxiesOn.next((Iterator *)&nextP));
      }
      
      nextC.id++;
      c = (computeInfo *) donor->computeSet.next((Iterator *)&nextC);
    } // end of compute loop

#define REASSIGN(GRID) if (GRID->c) { deAssign(GRID->c, donor); \
        assign(GRID->c, GRID->p); bestP = GRID->p; }

    bestP = 0;
    // see if we have found a compute processor pair
    REASSIGN(bestPe[3])
    else REASSIGN(bestPe[4])
    else REASSIGN(bestPe[5])
#if USE_TOPOMAP
    else REASSIGN(goodPe[3])
    else REASSIGN(goodPe[4])
    else REASSIGN(goodPe[5])
#endif
    else REASSIGN(bestPe[1])
    else REASSIGN(bestPe[2])
    else REASSIGN(bestPe[0])
#if USE_TOPOMAP
    else REASSIGN(goodPe[1])
    else REASSIGN(goodPe[2])
    else REASSIGN(goodPe[0])
#endif

    if(bestP) {
      if(bestP->load > averageLoad) lightPes->remove(bestP);
      if(donor->load > overLoad*averageLoad)
        heavyPes->insert((InfoRecord *) donor);
      else
	lightPes->insert((InfoRecord *) donor);
      //CkPrintf("1st try succeeded\n");
      
      continue;
    }
    //else
      //CkPrintf("1st try failed\n");

#if USE_TOPOMAP
    // if this fails, look at the inner brick
    int found = 0;
    int p1, p2, pe, x1, x2, xm, xM, y1, y2, ym, yM, z1, z2, zm, zM, t1, t2;
    int dimNX, dimNY, dimNZ, dimNT;
    double minLoad;

    good.c = 0; good.p = 0;
    minLoad = overLoad*averageLoad;
    nextC.id = 0;
    c = (computeInfo *)donor->computeSet.iterator((Iterator *)&nextC);
 
    while(c) {    
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

      // to expand the inner brick by some hops
#if 0
      if(xm>=EXPAND_INNER_BRICK) xm=xm-EXPAND_INNER_BRICK; else xm=0;
      if(ym>=EXPAND_INNER_BRICK) ym=ym-EXPAND_INNER_BRICK; else ym=0;
      if(zm>=EXPAND_INNER_BRICK) zm=zm-EXPAND_INNER_BRICK; else zm=0;

      xM=xM+EXPAND_INNER_BRICK;
      yM=yM+EXPAND_INNER_BRICK;
      zM=zM+EXPAND_INNER_BRICK;
#endif

      // first go over the processors inside the brick and choose the least 
      for(int i=xm; i<=xM; i++)
        for(int j=ym; j<=yM; j++)
	  for(int k=zm; k<=zM; k++)
	    for(int l=0; l<dimNT; l++)
	    {
	      pe = tmgr.coordinatesToRank(i%dimNX, j%dimNY, k%dimNZ, l);
	      p = &processors[pe];
	      if(c->load + p->load < minLoad) {
                minLoad = c->load + p->load;
	        good.c = c;
	        good.p = p;
	      }
	    }
      nextC.id++;
      c = (computeInfo *) donor->computeSet.next((Iterator *)&nextC);
    }

    if(good.c) {
      found = 1;
      //CkPrintf("2nd try succeeded\n");
    }
    else {
      found = 0;
      //CkPrintf("2nd try failed\n");
    }

    // if that also fails, look at the outer brick
    minLoad = overLoad * averageLoad;
    if(found==0) {
      good.c = 0; good.p = 0;
      p = 0;

      nextC.id = 0;
      c = (computeInfo *)donor->computeSet.iterator((Iterator *)&nextC);
      while(c) {
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

        for(int i=xM+1; i<xm+dimNX; i++)
          for(int j=0; j<dimNY; j++)
	    for(int k=0; k<dimNZ; k++)
	      for(int l=0; l<dimNT; l++)
	      {
	        pe = tmgr.coordinatesToRank(i%dimNX, j%dimNY, k%dimNZ, l);
	        p = &processors[pe];
	        if(c->load + p->load < minLoad) {
	          good.c = c;
	          good.p = p;
		  found = 1; break;
	        }
	      }

	if(found==1)
	  break;
	else {
          for(int j=yM+1; j<ym+dimNY; j++)
	    for(int i=xm; i<=xM; i++)
	      for(int k=0; k<dimNZ; k++)
	        for(int l=0; l<dimNT; l++)
	        {
		  pe = tmgr.coordinatesToRank(i%dimNX, j%dimNY, k%dimNZ, l);
		  p = &processors[pe];
		  if(c->load + p->load < minLoad) {
		    good.c = c;
		    good.p = p;
		    found = 1; break;
		  }
	        }
	}
 
	if(found==1)
	  break;
	else {
	  for(int k=zM+1; k<zm+dimNZ; k++)
	    for(int i=xm; i<=xM; i++)
	      for(int j=ym; j<=yM; j++)
	        for(int l=0; l<dimNT; l++)
	        {
		  pe = tmgr.coordinatesToRank(i%dimNX, j%dimNY, k%dimNZ);
		  p = &processors[pe];
		  if(c->load + p->load < minLoad) {
		    good.c = c;
		    good.p = p;
		    found = 1; break;
		  }
	        }
	}

	if(found==1) break;

	nextC.id++;
	c = (computeInfo *) donor->computeSet.next((Iterator *)&nextC);
      } 
    }

    if(found == 1) {
      deAssign(good.c, donor);
      assign(good.c, good.p);
      if (good.p->load > averageLoad) lightPes->remove(good.p);
      if (donor->load > overLoad*averageLoad)
        heavyPes->insert((InfoRecord *) donor);
      else
        lightPes->insert((InfoRecord *) donor);
      continue;
    }
    else {
      done = 0;
      break;
    }

#else

    int found=0;
    // find the first processor to place the compute on
    p = (processorInfo *)lightPes->iterator((Iterator *) &nextP);
    if(found == 0) {
     while (p)
      {
        nextC.id = 0;
        c = (computeInfo *)donor->computeSet.iterator((Iterator *)&nextC);
        while (c)
        {
          /*if(c->load + p->load < overLoad*averageLoad)
          {
	    good.c = c;
	    good.p = p;
            found = 1;
            break;
          }*/
	  selectPes(p, c);
          nextC.id++;
          c = (computeInfo *) donor->computeSet.next((Iterator *)&nextC);
        }
        //if(found == 1)
        //  break;
        p = (processorInfo *)lightPes->next((Iterator *) &nextP);
      }

      bestP = 0;
      REASSIGN(bestPe[3])
      else REASSIGN(bestPe[4])
      else REASSIGN(bestPe[5])
      else REASSIGN(bestPe[1])
      else REASSIGN(bestPe[2])
      else REASSIGN(bestPe[0])
    }

    if(bestP) {
      if(bestP->load > averageLoad) lightPes->remove(bestP);
      if(donor->load > overLoad*averageLoad)
        heavyPes->insert((InfoRecord *) donor);
      else
	lightPes->insert((InfoRecord *) donor);
      continue;
    }
    else {
      done = 0;
      break;
    }
#endif // USE_TOPOMAP
 
  } // end of while loop

  delete heavyPes;
  delete lightPes;

  return done;
}

void RefineTorusLB::selectPes(processorInfo *p, computeInfo *c) {
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
      pcpair* &newp = bestPe[index];

      if (!(newp->c) || ((p->load + c->load) < (newp->p->load + newp->c->load))) {
        newp->p = p;
        newp->c = c;
      } 
#if USE_TOPOMAP
    } else {
      pcpair* &newp = goodPe[index];

      if (!(newp->c) || ((p->load + c->load) < (newp->p->load + newp->c->load))) {
        newp->p = p;
        newp->c = c;

      }
    }
#endif
  }
}

