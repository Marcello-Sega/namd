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
  processorInfo *bestP, *bestP0, *bestP1, *bestP2;
  bestP0 = bestP1 = bestP2 = (processorInfo *)0;
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
  int i;
  overLoad = 1.2;
  for (i=0; i<numComputes; i++) {
    c = (computeInfo *) computesHeap->deleteMax();
    if (c->processor != -1) continue; // skip to the next compute;
    heapIterator nextProcessor;
    processorInfo *p = (processorInfo *) 
      pes->iterator((heapIterator *) &nextProcessor);
    // bestSize0 = bestSize1 = bestSize2 = 0;
    bestP0 = bestP1 = bestP2 = (processorInfo *)0;
    while (p) {
      int n=0;
      n = numAvailable(c,p);
      switch(n){
      case 0:
	if (!bestP0) {
	  if (c->load + p->load < overLoad*averageLoad)
	    bestP0 = p;
	} else {
	  if ( ( c->load + p->load < overLoad*averageLoad) && 
	       (p->load<bestP0->load))
	    bestP0 = p;
	}
	break;
      case 1:
	if (!bestP1) {
	  if (c->load + p->load < overLoad*averageLoad)
	    bestP1 = p;
	} else {
	  if (( c->load + p->load < overLoad*averageLoad) &&
	      (p->load < bestP1->load))
	    bestP1 = p;
	}
	break;
      case 2: 
	if (!bestP2) {
	  if (c->load + p->load < overLoad*averageLoad)
	    bestP2 = p;
	} else {
	  if (( c->load + p->load < overLoad*averageLoad) &&
	      (p->load < bestP2->load))
	    bestP2 = p;
	}
	break;
      default:
	iout << iINFO  << "Error. Illegal number of proxies.\n" << endi;    
      }
      p = (processorInfo *) pes->next(&nextProcessor);
    }

    if (numAssigned >= 0) {
      if (bestP2) {
	assign(c, bestP2);
	numAssigned++;
	numAssignedP2++;
      } else if (bestP1) {
	assign(c, bestP1);
	numAssigned++;
	numAssignedP1++;
      } else if (bestP0) {
        assign(c, bestP0);
        numAssigned++;
        numAssignedP0++;
      } else { 
        iout << iINFO  << "Alg 7 No receiver found 1" << "\n" <<endi;
        break;
      }
    } else {
      // At start, load is most important, rather than communications
      int *numAssignedptr = &numAssignedP2;
      bestP = bestP2;
      if (!bestP || (bestP1 && (bestP1->load < 0.8 * bestP->load)) ) {
	bestP=bestP1;
	numAssignedptr = &numAssignedP1;
      }
      if (!bestP || (bestP0 && (bestP0->load < 0.75 * bestP->load))) {
	bestP=bestP0;
	numAssignedptr = &numAssignedP0;
      }
      if (!bestP) {
	iout << iINFO  << "Alg7 No receiver found 2" << "\n" <<endi;
	break;
      }
      assign(c,bestP);
      (*numAssignedptr)++;
      numAssigned++;
      numAssignedP4++;
    }
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















