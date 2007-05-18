/**
***  Copyright (c) 1995, 1996, 1997, 1998, 1999, 2000 by
***  The Board of Trustees of the University of Illinois.
***  All rights reserved.
**/

#ifndef REFINEONLY_DEFS_H
#define REFINEONLY_DEFS_H

#include "elements.h"
#include "heap.h"

// #define LDB_DEBUG

class Rebalancer {
private:
  int bytesPerAtom;
  void InitProxyUsage();

  struct pcpair {
    processorInfo *p;
    computeInfo *c;
    pcpair() : p(0),c(0) {;}
  };
  typedef pcpair pcgrid[3][3][2];

  void refine_togrid(pcgrid &grid, double thresholdLoad,
                        processorInfo *p, computeInfo *c);

protected: 
  const char *strategyName;
  computeInfo *computes;
  patchInfo *patches;
  processorInfo *processors;
  minHeap *pes;
  maxHeap *computePairHeap;
  maxHeap *computeSelfHeap;
  maxHeap *computeBgPairHeap;
  maxHeap *computeBgSelfHeap;
  int P;
  int numPatches;
  int numComputes;
  int numProxies;
  int numPesAvailable;
  double averageLoad;
  int isAvailableOn(patchInfo *patch, processorInfo *p);
  void numAvailable(computeInfo *c, processorInfo *p,
           int *nPatches, int *nProxies, int *isBadForCommunication);

  void strategy();
  void makeHeaps();
  void assign(computeInfo *c, processorInfo *pRec);
  void assign(computeInfo *c, int p);
  void deAssign(computeInfo *c, processorInfo *pRec);
  int refine();
  void multirefine(double overload_start=1.02);
  void printSummary();
  void printResults();
  void printLoads();
  double computeAverage();
  void adjustBackgroundLoadAndComputeAverage();
  double computeMax();
  double overLoad;

public:
  Rebalancer(computeInfo *computeArray, patchInfo *patchArray,
             processorInfo *processorArray,
             int nComps, int nPatches, int nPes);
  ~Rebalancer();
};

#if CHARM_VERSION > 50913 && USE_TOPOMAP 
#include "TopoManager.h"
#endif
#endif
