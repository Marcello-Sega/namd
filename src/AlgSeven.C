/**
***  Copyright (c) 1995, 1996, 1997, 1998, 1999, 2000 by
***  The Board of Trustees of the University of Illinois.
***  All rights reserved.
**/

#include <iostream.h>
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

  int i,j;
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
  for (i=0; i<numComputes; i++) {
    c = (computeInfo *) computesHeap->deleteMax();
    if (c->processor != -1) continue; // skip to the next compute;
    heapIterator nextProcessor;
    processorInfo *p = (processorInfo *) 
      pes->iterator((heapIterator *) &nextProcessor);
    for(i=0;i<3;i++)
      for(j=0;j<3;j++)
	goodP[i][j]=0;
    while (p) {
      int nPatches = numPatchesAvail(c,p);
      int nProxies = numProxiesAvail(c,p);
      if (nPatches < 0 || nPatches > 2)
	iout << iERROR << "Too many patches: " << nPatches << "\n" << endi;
      if (nProxies < 0 || nProxies > 2)
	iout << iERROR << "Too many proxies: " << nProxies << "\n" << endi;

      if (!goodP[nPatches][nProxies]) {
	if (nPatches == 0 && nProxies == 0)
	  goodP[0][0] = p;
	else if (c->load + p->load < overLoad*averageLoad)
	  goodP[nPatches][nProxies] = p;
      } else {
	if (( c->load + p->load < overLoad*averageLoad) &&
	    (p->load < goodP[nPatches][nProxies]->load))
	  goodP[nPatches][nProxies] = p;
      }
      p = (processorInfo *) pes->next(&nextProcessor);
    }

    //    if (numAssigned >= 0) {  Else is commented out below

    processorInfo* selectedP;
    if (goodP[2][0]) {
      // Two home, no proxies
      assign(c, goodP[2][0]);
      numAssigned++;
    } else if (goodP[1][1]) {
      // One home, one proxy
      assign(c, goodP[1][1]);
      numAssigned++;
    } else if (goodP[0][2]) {
      // No home, two proxies
      assign(c, goodP[0][2]);
      numAssigned++;
    } else if (goodP[1][0]) {
      // One home, no proxies
      assign(c, goodP[1][0]);
      numAssigned++;
    } else if (goodP[0][1]) {
      // No home, one proxy
      assign(c, goodP[0][1]);
      numAssigned++;
    } else if (goodP[0][0]) {
      // No home, no proxies
      assign(c, goodP[0][0]);
      numAssigned++;
    } else {
      iout << iERROR  << "*** Alg 7 No receiver found 1 ***" << "\n" <<endi;
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


  // The New refinement procedure.  This is identical to the code in
  // RefineOnly.C, and probably should be merged with that code to form
  // a binary-search function

  double avg = computeAverage();
  double max = computeMax();

  const double overloadStep = 0.01;
  const double overloadStart = 1.02;
  double dCurOverload = max / avg;

  int minOverload = 0;
  int maxOverload = (int)((dCurOverload - overloadStart)/overloadStep + 1);
  double dMinOverload = minOverload * overloadStep + overloadStart;
  double dMaxOverload = maxOverload * overloadStep + overloadStart;

  iout << iINFO
       << "Balancing from " << minOverload << " = " << dMinOverload 
       << " to " << maxOverload << "=" << dMaxOverload 
       << " dCurOverload=" << dCurOverload << " max=" << max << " avg=" << avg
       << "\n" << endi;

  int curOverload;
  int refineDone = 0;

  overLoad = dMinOverload;
  if (refine())
    refineDone = 1;
  else {
    overLoad = dMaxOverload;
    if (!refine()) {
      iout << iINFO << "ERROR: Could not refine at max overload\n" << endi;
      refineDone = 1;
    }
  }

  // Scan up, until we find a refine that works
  while (!refineDone) {
    if (maxOverload - minOverload <= 1)
      refineDone = 1;
    else {
      curOverload = (maxOverload + minOverload ) / 2;

      overLoad = curOverload * overloadStep + overloadStart;
      iout << iINFO << "Testing curOverload " << curOverload 
	   << "=" << overLoad << " [min,max]=" 
	   << minOverload << ", " << maxOverload
	   << "\n" << endi;
      if (refine())
	maxOverload = curOverload;
      else
	minOverload = curOverload;
    }
  }

  //  iout << iINFO  << "num assigned: " << numAssigned << endi;
  //  iout << iINFO  << "Starting overLoad = " << overLoad << endi;
  //  iout << iINFO  << "Ending overLoad = " << overLoad << endi;
  //  iout << iINFO
  //   << "After assignment\n" << endi;
  // printLoads();
}















