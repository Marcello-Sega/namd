/**
***  Copyright (c) 1995, 1996, 1997, 1998, 1999, 2000 by
***  The Board of Trustees of the University of Illinois.
***  All rights reserved.
**/

#ifndef REFINEONLY_DEFS_H
#define REFINEONLY_DEFS_H

#include "elements.h"
#include "heap.h"


class Rebalancer {
private:
  int bytesPerAtom;
  void InitProxyUsage();

protected: 
  const char *strategyName;
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
  int numPatchesAvail(computeInfo *c, processorInfo *p);
  int numProxiesAvail(computeInfo *c, processorInfo *p);

  void strategy();
  void makeHeaps();
  void assign(computeInfo *c, processorInfo *pRec);
  void assign(computeInfo *c, int p);
  void deAssign(computeInfo *c, processorInfo *pRec);
  int refine();
  void multirefine();
  void printResults();
  void printLoads();
  double computeAverage();
  double computeMax();
  double overLoad;

public:
  Rebalancer(computeInfo *computeArray, patchInfo *patchArray,
             processorInfo *processorArray,
             int nComps, int nPatches, int nPes);
  ~Rebalancer();
};


#endif
