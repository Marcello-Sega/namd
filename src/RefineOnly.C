/**
***  Copyright (c) 1995, 1996, 1997, 1998, 1999, 2000 by
***  The Board of Trustees of the University of Illinois.
***  All rights reserved.
**/

#include "InfoStream.h"
#include "RefineOnly.h"

RefineOnly::RefineOnly(computeInfo *computeArray, patchInfo *patchArray, 
	   processorInfo *processorArray, int nComps, 
	   int nPatches, int nPes) :
Rebalancer(computeArray, patchArray, 
	   processorArray, nComps, 
	   nPatches, nPes)
{
strategyName = "RefineOnly";
strategy();
}

void RefineOnly::strategy()
{ 
iout << iINFO << "numComputes: " << numComputes << "\n" << endi;
  for (int i=0; i<numComputes; i++)
    assign((computeInfo *) &(computes[i]),
	   (processorInfo *) &(processors[computes[i].oldProcessor]));
	 
  // merge refinement code and move to Rebalancer
  multirefine();

#if 0
  double avg = computeAverage();
  double max = computeMax();

#if 0
  iout << iINFO
       << "------------------------------------------------------------\n"
       << iINFO << "Before load balancing (measured stats):\n" << endi;
#endif
  printLoads();

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
#endif
    
  //  while (!refine())
  //    overLoad += .01;

  computeAverage();
#if 0
  iout << iINFO
       << "------------------------------------------------------------\n"
       << iINFO 
       << "After load balancing (predicted stats):\n" << endi;
#endif
  printLoads();
#if 0
       << "------------------------------------------------------------\n"
       << endi;
#endif
  double max = computeMax();
}
