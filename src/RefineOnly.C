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
  iout << iINFO
       << "------------------------------------------------------------\n"
       << iINFO << "Before load balancing (measured stats):\n" << endi;
  printLoads();

  double max = computeMax();
  refine();

  computeAverage();
  iout << iINFO
       << "------------------------------------------------------------\n"
       << iINFO 
       << "After load balancing (predicted stats):\n" << endi;
  printLoads();
  iout << iINFO 
       << "------------------------------------------------------------\n"
       << endi;
  max = computeMax();
}
