/**
***  Copyright (c) 1995, 1996, 1997, 1998, 1999, 2000 by
***  The Board of Trustees of the University of Illinois.
***  All rights reserved.
**/

#ifndef COMPUTENONBONDEDSELF_H
#define COMPUTENONBONDEDSELF_H

#include "ComputePatch.h"
#include "ComputeNonbondedUtil.h"

class ComputeNonbondedSelf : public ComputePatch, private ComputeNonbondedUtil {

public:
  ComputeNonbondedSelf(ComputeID c, PatchID pid,
	int minPartition = 0, int maxPartition = 1, int numPartitions = 1);
  virtual ~ComputeNonbondedSelf();

protected :
  virtual void initialize();
  virtual void doForce(Position* p, Results* r, AtomProperties* a);

  PositionBox<Patch> *avgPositionBox;

  SubmitReduction *reduction;

  int minPart, maxPart, numParts;

};

#endif

