/**
***  Copyright (c) 1995, 1996, 1997, 1998, 1999, 2000 by
***  The Board of Trustees of the University of Illinois.
***  All rights reserved.
**/

#include <iostream.h>
#include "common.h"
#include "InfoStream.h"
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
  // int numAssigned1 = 0;
  // int numAssigned2 = 0;
  int numAssignedP2 = 0;
  int numAssignedP1 = 0;
  int numAssignedP0 = 0;
  int numAssignedP4 = 0;

  //   for (int i=0; i<numPatches; i++)
  //     { cout << "(" << patches[i].Id << "," << patches[i].processor ;}
  overLoad = 1.2;
  for (int ic=0; ic<numComputes; ic++) {
    c = (computeInfo *) computesHeap->deleteMax();
    if ( ! c ) NAMD_bug("Alg7: computesHeap empty!");
    if (c->processor != -1) continue; // skip to the next compute;
    heapIterator nextProcessor;
    processorInfo *p = (processorInfo *) 
      pes->iterator((heapIterator *) &nextProcessor);
    int i,j;
    for(i=0;i<3;i++)
      for(j=0;j<3;j++) {
	goodP[i][j]=0;
	poorP[i][j]=0;
      }
    while (p) {
      int nPatches = numPatchesAvail(c,p);
      int nProxies = numProxiesAvail(c,p);
      if (nPatches < 0 || nPatches > 2)
	iout << iERROR << "Too many patches: " << nPatches << "\n" << endi;
      if (nProxies < 0 || nProxies > 2)
	iout << iERROR << "Too many proxies: " << nProxies << "\n" << endi;

      if (!goodP[nPatches][nProxies] ||
	    (p->load < goodP[nPatches][nProxies]->load)) {
        if (c->load + p->load < overLoad*averageLoad) {
	  goodP[nPatches][nProxies] = p;
        }
      }
      if (!poorP[nPatches][nProxies] ||
	    (p->load < poorP[nPatches][nProxies]->load)) {
	poorP[nPatches][nProxies] = p;   // fallback
      }
      p = (processorInfo *) pes->next(&nextProcessor);
    }

    //    if (numAssigned >= 0) {  Else is commented out below

    p = 0;
    if ((p = goodP[2][0])    // Two home, no proxies
     || (p = goodP[1][1])    // One home, one proxy
     || (p = goodP[0][2])    // No home, two proxies
     || (p = goodP[1][0])    // One home, no proxies
     || (p = goodP[0][1])    // No home, one proxy
     || (p = goodP[0][0])    // No home, no proxies
     || (p = poorP[2][0])    // Two home, no proxies, overload
     || (p = poorP[1][1])    // One home, one proxy, overload
     || (p = poorP[0][2])    // No home, two proxies, overload
     || (p = poorP[1][0])    // One home, no proxies, overload
     || (p = poorP[0][1])    // No home, one proxy, overload
     || (p = poorP[0][0])    // No home, no proxies, overload
       ) {
      assign(c,p); numAssigned++;
    } else {
      NAMD_bug("*** Alg 7 No receiver found 1 ***");
      break;
    }

//     } else {
//       // At start, load is most important, rather than communications
//       int *numAssignedptr = &numAssignedP2;
//       bestP = bestP2;
//       if (!bestP || (bestP1 && (bestP1->load < 0.8 * bestP->load)) ) {
// 	bestP=bestP1;
// 	numAssignedptr = &numAssignedP1;
//       }
//       if (!bestP || (bestP0 && (bestP0->load < 0.75 * bestP->load))) {
// 	bestP=bestP0;
// 	numAssignedptr = &numAssignedP0;
//       }
//       if (!bestP) {
// 	iout << iERROR  << "*** Alg7 No receiver found 2 ***" << "\n" <<endi;
// 	break;
//       }
//       assign(c,bestP);
//       (*numAssignedptr)++;
//       numAssigned++;
//       numAssignedP4++;
//     }
  }

#ifdef DEBUG
  iout << iINFO
       << "numAssigned = " << numAssigned
       << "\nnumAssignedP2 = " << numAssignedP2
       << "\nnumAssignedP1 = " << numAssignedP1
       << "\nnumAssignedP0 = " << numAssignedP0
       << "\nnumAssignedP4 = " << numAssignedP4
       << "\n" << endi;
#endif

#if 0
  // The old Refinement procedure

  overLoad = 1.02;
  for (; !refine(); overLoad += .01);
#endif


  // binary-search refinement procedure
  multirefine();


  //  iout << iINFO  << "num assigned: " << numAssigned << endi;
  //  iout << iINFO  << "Starting overLoad = " << overLoad << endi;
  //  iout << iINFO  << "Ending overLoad = " << overLoad << endi;
  //  iout << iINFO
  //   << "After assignment\n" << endi;
  // printLoads();

  CmiPrintf("Alg7 finish time: %f.\n", CmiWallTimer()-startTime);
}















