#include <iostream.h>
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
  cout << "Average Load is " << averageLoad << "\n";
  cout << "Maximum Load is " << max << "\n";
  refine();
  computeAverage();
  printLoads();
  max = computeMax();
  cout << "Average Load is " << averageLoad << "\n";
  cout << "Maximum Load is " << max << "\n";
}
