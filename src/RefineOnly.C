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
  for (int i=0; i<numComputes; i++)
    assign((computeInfo *) &(computes[i]),
	   (processorInfo *) &(processors[computes[i].oldProcessor]));
	 

  computeAverage();
#if 0
  iout << iINFO
       << "------------------------------------------------------------\n"
       << iINFO << "Before load balancing (measured stats):\n" << endi;
  printLoads();
#endif

  double max = computeMax();

  overLoad = 1.02;
  while (!refine())
    overLoad += .01;

  computeAverage();
#if 0
  iout << iINFO
       << "------------------------------------------------------------\n"
       << iINFO 
       << "After load balancing (predicted stats):\n" << endi;
  printLoads();
  iout << iINFO 
       << "------------------------------------------------------------\n"
       << endi;
#endif
  max = computeMax();
}
