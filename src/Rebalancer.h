/**
***  Copyright (c) 1995, 1996, 1997, 1998, 1999, 2000 by
***  The Board of Trustees of the University of Illinois.
***  All rights reserved.
**/

#ifndef REFINEONLY_DEFS_H
#define REFINEONLY_DEFS_H

#include "elements.h"
#include "heap.h"
#if CHARM_VERSION > 50913 && USE_TOPOMAP 
#include "TopoManager.h"
#endif
#include "ProxyMgr.decl.h"
#include "ProxyMgr.h"

// #define LDB_DEBUG

#include "ckhashtable.h"

#if USE_TOPOMAP
#include "TopoManager.h"
#endif

class ProxyUsageKey {
 protected:
  int      processor;
  int      patch;
  
 public:
  ProxyUsageKey (int pe, int patch) {
    this->processor = pe;
    this->patch     = patch;
  }

  CkHashCode hash (void) const {
    return (patch << 16) + processor;
  }

  static CkHashCode  staticHash (const void *x, size_t size) {
    return ((ProxyUsageKey *)x)->hash();
  }

  int compare (const ProxyUsageKey &in) const {
    if ((in.patch == patch) && (in.processor == processor))
      return 1;
    
    return 0;
  }
   
  static int staticCompare (const void *a, const void *b, size_t size) {
    return ((ProxyUsageKey *)a)->compare(* (ProxyUsageKey *)b);
  }
};

class ProxyUsage {
 protected:
  CkHashtableT <ProxyUsageKey, int>  htable;
  
 public:
  
  ProxyUsage () : htable (1217, 0.5) {}   //pass in a large prime close to 
                                          //1k processors

  void increment (int pe, int patch) {
    ProxyUsageKey  pkey (pe, patch);

    int val = htable.get (pkey);
    htable.put (pkey) =  val + 1;      
  }

  void decrement (int pe, int patch) {
    ProxyUsageKey  pkey (pe, patch);
    
    int val = htable.get (pkey);
    CkAssert (val > 0);
    val --;

    if (val == 0)
      htable.remove (pkey);
    else 
      htable.put (pkey) = val; 
  }
     
  int getVal (int pe, int patch) {
    ProxyUsageKey  pkey (pe, patch);  
    return htable.get (pkey);
  }
};
  

class Rebalancer {
public:
  struct pcpair {
    processorInfo *p;
    computeInfo *c;
    pcpair() : p(0),c(0) {;}
    void reset () { p = 0; c = 0; }
  };

private:
  int bytesPerAtom;
  
  typedef pcpair pcgrid[3][3][2];

  void refine_togrid(pcgrid &grid, double thresholdLoad,
                        processorInfo *p, computeInfo *c);

  ProxyUsage  proxyUsage;


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
  void makeTwoHeaps();
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
  void createSpanningTree();
  void decrSTLoad();
  void incrSTLoad();
  void InitProxyUsage();
#if USE_TOPOMAP
  TopoManager tmgr;
#endif

public:
  Rebalancer(computeInfo *computeArray, patchInfo *patchArray,
             processorInfo *processorArray,
             int nComps, int nPatches, int nPes);
  ~Rebalancer();
};

#endif
