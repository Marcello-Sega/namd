/**
***  Copyright (c) 1995, 1996, 1997, 1998, 1999, 2000 by
***  The Board of Trustees of the University of Illinois.
***  All rights reserved.
**/

#include <iostream.h>
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
      if (nPatches < 0 || nPatches > 2)
	iout << iERROR << "Too many patches: " << nPatches << "\n" << endi;
      if (nProxies < 0 || nProxies > 2)
	iout << iERROR << "Too many proxies: " << nProxies << "\n" << endi;

      if (c->load + p->load < overLoad*averageLoad) {
        processorInfo* &altp = goodP[nPatches][nProxies];
        if ( nPatches + nProxies == 2 ) {
          if (!altp || p->load < altp->load ) {
	    altp = p;
          }
        } else if (p->proxies.numElements() <
		((double)numProxies / (double)numPesAvailable + 2) ) {
          if (!altp || p->proxies.numElements() <
			altp->proxies.numElements() ) {
	    altp = p;
          }
        }
      }

      {
        processorInfo* &altp = poorP[nPatches][nProxies];
        if ( nPatches + nProxies == 2 ) {
          if (!altp || p->load < altp->load ) {
	    altp = p;
          }
        } else if (p->proxies.numElements() <
		((double)numProxies / (double)numPesAvailable + 2) ) {
          if (!altp || p->proxies.numElements() <
			altp->proxies.numElements() ) {
	    altp = p;
          }
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
     printLoads();
	      
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

  CmiPrintf("Alg7 finish time: %f.\n", CmiWallTimer()-startTime);
}















