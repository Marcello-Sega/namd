#ifndef REFINEONLY_DEFS_H
#define REFINEONLY_DEFS_H

class minheap;
class maxheap;

#include "elements.h"
#include "heap.h"


class Rebalancer 
{
protected: 
  int bytesPerAtom;
  computeInfo *computes;
  patchInfo *patches;
  processorInfo *processors;
  minHeap *pes;
  maxHeap *computesHeap;
  int P;
  int numPatches;
  int numComputes;
  double averageLoad;
  int isAvailableOn(patchInfo *patch, processorInfo *p);
  int numAvailable(computeInfo *c, processorInfo *p);

  void strategy();
  void makeHeaps();
  void assign(computeInfo *c, processorInfo *pRec);
  void assign(computeInfo *c, int p);
  void deAssign(computeInfo *c, processorInfo *pRec);
  int refine();
  void printResults();
  void printLoads();
  void computeAverage();
  double computeMax();
  
  void InitProxyUsage();


public:
  double overLoad;
  char *strategyName;
  Rebalancer() {}
  ~Rebalancer();
  Rebalancer(computeInfo *computeArray, patchInfo *patchArray,
             processorInfo *processorArray,
             int nComps, int nPatches, int nPes);
};


#endif
