/**
***  Copyright (c) 1995, 1996, 1997, 1998, 1999, 2000 by
***  The Board of Trustees of the University of Illinois.
***  All rights reserved.
**/

#ifndef ALGROB_H
#define ALGROB_H

//#include "elements.h"
#include "PatchMap.h"
#include "Rebalancer.h"

class AlgRecBisection : public Rebalancer 
{
private: 

typedef struct {
  int refno;
  double load;				// total load in this set
  int origin[3];
  int corner[3];
  int  count;
} Partition;

typedef struct {
  int id;
  int v[3];
  double load;
  int  refno;
  int  tv;
} ComputeLoad;

public:

typedef struct {
  int v;
  int id;
} VecArray;

private:
enum {XDIR=0, YDIR, ZDIR};

ComputeLoad *computeLoad;
VecArray  *(vArray[3]);
Partition *partitions;
Partition top_partition;
int npartition;
int currentp, refno;

void strategy();
void rec_divide(int, Partition&);
void setVal(int x, int y, int z);
int sort_partition(int x, int p, int r);
void qsort(int x, int p, int r);
void quicksort(int x);


public:
AlgRecBisection(computeInfo *computeArray, patchInfo *patchArray, 
	   processorInfo *processorArray, int nComps, 
	   int nPatches, int nPes);
};

#endif




