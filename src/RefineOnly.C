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
{ printLoads();
   for (int i=0; i<numComputes; i++)
    assign((computeInfo *) &(computes[i]),
	   (processorInfo *) &(processors[computes[i].oldProcessor]));
	 

  computeAverage();
  printLoads();
  double max = computeMax();
  iout << "Average Load is " << averageLoad << "\n";
  iout << "Maximum Load is " << max << "\n" << endi;
  refine();
  computeAverage();
  printLoads();
  max = computeMax();
  iout << "Average Load is " << averageLoad << "\n";
  iout << "Maximum Load is " << max << "\n" << endi;
}
