#ifndef ALG7_H
#define ALG7_H

//#include "elements.h"
#include "Rebalancer.h"

class Alg7 : public Rebalancer 
{
private: 
void strategy();


public:
Alg7(computeInfo *computeArray, patchInfo *patchArray, 
	   processorInfo *processorArray, int nComps, 
	   int nPatches, int nPes);
};

#endif




