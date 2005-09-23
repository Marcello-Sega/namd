/**
***  Copyright (c) 1995, 1996, 1997, 1998, 1999, 2000 by
***  The Board of Trustees of the University of Illinois.
***  All rights reserved.
**/

#include "common.h"
#include "InfoStream.h"
#include "Node.h"
#include "Alg7.h"

#define TINYLOAD 0.0005

Alg7::Alg7(computeInfo *computeArray, patchInfo *patchArray, 
	   processorInfo *processorArray, int nComps, 
	   int nPatches, int nPes) :
  Rebalancer(computeArray, patchArray, 
	     processorArray, nComps, 
	     nPatches, nPes)
{
strategyName = "Alg7";
strategy();
}

extern int isPmeProcessor(int);

void Alg7::togrid(processorInfo* goodP[3][3], processorInfo* poorP[3][3],
			processorInfo *p, computeInfo *c) {
      if(p->available == CmiFalse) return;

      int nPatches = numPatchesAvail(c,p);
      int nProxies = numProxiesAvail(c,p);
      if ( nProxies < 0 ) { nPatches = nProxies = 0; }
      if ( nPatches + nProxies < 2 && p->proxies.numElements() > 6 &&
          p->proxies.numElements() >
		((double)numProxies / (double)numPesAvailable + 3) ) {
        nPatches = nProxies = 0;
      }

      if (nPatches < 0 || nPatches > 2)
	iout << iERROR << "Too many patches: " << nPatches << "\n" << endi;
      if (nProxies < 0 || nProxies > 2)
	iout << iERROR << "Too many proxies: " << nProxies << "\n" << endi;

      if (c->load + p->load < overLoad*averageLoad) {
        processorInfo* &altp = goodP[nPatches][nProxies];	

#if CMK_VERSION_BLUEGENE
	if(!altp)
	  altp = p;
	else {
	  //Find processors that are patch neighbors on the BGL torus
	  int neighbor = 0, neighbor_alt = 0;
	  
	  BGLTorusManager *tmgr = BGLTorusManager::getObject();
	  /*
	    if((tmgr->isNeighbor(altp->Id, patches[c->patch1].processor) ||
	    tmgr->isNeighbor(altp->Id, patches[c->patch2].processor)))
	    neighbor_alt = 1;
	    
	    if(tmgr->isNeighbor(p->Id, patches[c->patch1].processor) ||
	    tmgr->isNeighbor(p->Id, patches[c->patch2].processor))
	    neighbor = 1;
	  */
	  
	  if(tmgr->isNeighborOfBoth(altp->Id, patches[c->patch1].processor,
				    patches[c->patch2].processor, 1))
	    neighbor_alt = 1;
	  
	  if(tmgr->isNeighborOfBoth(p->Id, patches[c->patch1].processor, 
				    patches[c->patch2].processor, 1))
	    neighbor = 1;
	  
	  if(neighbor_alt == 1 && neighbor == 1) {
	    //Both are neighbors, only replace if lighter
	    if (p->load < altp->load ) {
	      altp = p;
	    }
	  }
	  else if(neighbor_alt == 0 && neighbor == 1)
	    //Previous was not a neighbor, kick him out
	    altp = p;
	  else if(neighbor_alt == 1 && neighbor == 0)
	    ;      //Give preference to good neighbors
	  else {
	    //Both not neighbors, choose nearby node to minimize hop bytes
	    /*
	      if (!altp || p->load < altp->load ) {
	      altp = p;
	      }
	    */

	    int alt_dist = 0, dist = 0;	    
	    int ax,ay,az, x,y,z, p1x,p1y,p1z, p2x,p2y,p2z;
	    
	    tmgr->getCoordinatesByRank(altp->Id, ax,ay,az);
	    tmgr->getCoordinatesByRank(p->Id, x,y,z);
	    
	    tmgr->getCoordinatesByRank(patches[c->patch1].processor, p1x,p1y,p1z);
	    tmgr->getCoordinatesByRank(patches[c->patch2].processor, p2x,p2y,p2z);
	    
	    alt_dist = abs(p1x - ax) + abs(p2x - ax) +
	      abs(p1y - ay) + abs(p1z - az) +
	      abs(p2y - ay) + abs(p2z - az);
	    
	    dist = abs(p1x - x) + abs(p2x - x) +
	      abs(p1y - y) + abs(p1z - z) +
	      abs(p2y - y) + abs(p2z - z);
	    
	    if(alt_dist > dist)
	      altp = p;	  
	  }
	}
#else 
        if (!altp || p->load < altp->load ) {	
	  altp = p;
        }
#endif	  
      }
      
      {
        processorInfo* &altp = poorP[nPatches][nProxies];
        if (!altp || p->load < altp->load ) {
	  altp = p;
        }
      }
}

void Alg7::strategy()
{
  // double bestSize0, bestSize1, bestSize2;
  computeInfo *c;
  int numAssigned;
  processorInfo* goodP[3][3];  // goodP[# of real patches][# of proxies]
  processorInfo* poorP[3][3];  // fallback option

  double startTime = CmiWallTimer();

  //   iout << iINFO  << "calling makeHeaps. \n";
  makeHeaps();
  computeAverage();
  //   iout << iINFO
  //	<< "Before assignment\n" << endi;
  //   printLoads();
	      
  numAssigned = 0;

  //   for (int i=0; i<numPatches; i++)
  //     { cout << "(" << patches[i].Id << "," << patches[i].processor ;}
  overLoad = 1.2;
  for (int ic=0; ic<numComputes; ic++) {
    c = (computeInfo *) computesHeap->deleteMax();
    if ( ! c ) NAMD_bug("Alg7: computesHeap empty!");
    if (c->processor != -1) continue; // skip to the next compute;
    int i,j;
    for(i=0;i<3;i++)
      for(j=0;j<3;j++) {
	goodP[i][j]=0;
	poorP[i][j]=0;
      }

    // first try for at least one proxy
    {
      Iterator nextProc;
      processorInfo *p;

      p = &processors[patches[c->patch1].processor];
      togrid(goodP, poorP, p, c);

      p = &processors[patches[c->patch2].processor];
      togrid(goodP, poorP, p, c);

      p = (processorInfo *)patches[c->patch1].
                            proxiesOn.iterator((Iterator *)&nextProc);
      while (p) {
        togrid(goodP, poorP, p, c);
        p = (processorInfo *)patches[c->patch1].
                            proxiesOn.next((Iterator*)&nextProc);
      }

      p = (processorInfo *)patches[c->patch2].
                            proxiesOn.iterator((Iterator *)&nextProc);
      while (p) {
        togrid(goodP, poorP, p, c);
        p = (processorInfo *)patches[c->patch2].
                            proxiesOn.next((Iterator*)&nextProc);
      }
      p = 0;
      if ((p = goodP[1][1])    // One home, one proxy
       || (p = goodP[2][0])    // Two home, no proxies
       || (p = goodP[0][2])    // No home, two proxies
       || (p = goodP[1][0])    // One home, no proxies
       || (p = goodP[0][1])    // No home, one proxy
       || (p = goodP[0][0])    // No home, no proxies
         ) {
        assign(c,p); numAssigned++;
        continue;
      }
    }

    // no luck, do it the long way

    heapIterator nextProcessor;
    processorInfo *p = (processorInfo *) 
      pes->iterator((heapIterator *) &nextProcessor);
    while (p) {
      togrid(goodP, poorP, p, c);
      p = (processorInfo *) pes->next(&nextProcessor);
    }

    //    if (numAssigned >= 0) {  Else is commented out below

    p = 0;
    if ((p = goodP[1][1])    // One home, one proxy
     || (p = goodP[2][0])    // Two home, no proxies
     || (p = goodP[0][2])    // No home, two proxies
     || (p = goodP[1][0])    // One home, no proxies
     || (p = goodP[0][1])    // No home, one proxy
     || (p = goodP[0][0])    // No home, no proxies
       ) {
      assign(c,p); numAssigned++;
   } else if (
        (p = poorP[2][0])    // Two home, no proxies, overload
     || (p = poorP[1][1])    // One home, one proxy, overload
     || (p = poorP[0][2])    // No home, two proxies, overload
     || (p = poorP[1][0])    // One home, no proxies, overload
     || (p = poorP[0][1])    // No home, one proxy, overload
     || (p = poorP[0][0])    // No home, no proxies, overload
       ) {
      iout << iWARN << "overload assign to " << p->Id << "\n" << endi;
      assign(c,p); numAssigned++;
    } else {
      NAMD_bug("*** Alg 7 No receiver found 1 ***");
      break;
    }

  }

  printLoads();

  // binary-search refinement procedure
  multirefine();
  printLoads();

  // CmiPrintf("Alg7 finish time: %f.\n", CmiWallTimer()-startTime);
}

