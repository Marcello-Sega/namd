#ifndef ELEMENTS_DEFS_H
#define ELEMENTS_DEFS_H

#include "Set.h"
class minHeap;
class maxHeap;

class InfoRecord
{
public:
  double load;
  int Id; // should replace other Ids.
};


class computeInfo : public InfoRecord
{
public: 
  /*   int computeId; replaced by Id */
  int patch1, patch2;
  int processor; // caller to ReBalancer MAY leave this field -1, 
  int oldProcessor; // stores the current assignment of the compute object.
};

class patchInfo : public InfoRecord
{
public:
  int processor;
  int numAtoms;
  //  int patchId; replaced by inherited "Id"
//  Set *computeSet; // caller to ReBalancer should leave this field NULL.
  Set *proxiesOn;  // caller to ReBalancer should fill in the forced proxies
};

class processorInfo: public InfoRecord
{
public:
//   int processorNum; replaced by inherited "Id".
  double backgroundLoad; // background work pre-assigned to the processor.
  double computeLoad;  //load due to computes. The total load is computed
			   // by adding these two.		     
  Set *patchSet; // caller to ReBalancer should leave this field NULL.
  Set *proxies;  // caller to ReBalancer should fill in the forced proxies
  Set *computeSet; // caller to ReBalancer should leave this field NULL.
  maxHeap *computesWithBoth; // heap of compute objects with both patches 
				  // (either in the home form, or proxy) here.
  maxHeap *computesWithOne;  // with one of the two here.
};


#endif
