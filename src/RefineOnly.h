#ifndef _REFINEONLY_H_
#define _REFINEONLY_H_

#include "elements.h"
#include "Rebalancer.h"

class RefineOnly : public Rebalancer 
{
private: 
void strategy();


public:
RefineOnly(computeInfo *computeArray, patchInfo *patchArray, 
	   processorInfo *processorArray, int nComps, 
	   int nPatches, int nPes);
};

#endif
