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
	ComputeNonbondedWorkArrays* _workArrays,
	int minPartition = 0, int maxPartition = 1, int numPartitions = 1);
  virtual ~ComputeNonbondedSelf();

protected :
  virtual void initialize();
  virtual int noWork();
  virtual void doForce(CompAtom* p, CompAtomExt* pExt, Results* r);
  Box<Patch,CompAtom> *avgPositionBox;
  // BEGIN LA
  Box<Patch,CompAtom> *velocityBox;
  // END LA

  SubmitReduction *reduction;
  SubmitReduction *pressureProfileReduction;
  BigReal *pressureProfileData;

  ComputeNonbondedWorkArrays* const workArrays;

  Pairlists pairlists;
  int pairlistsValid;
  BigReal pairlistTolerance;

  int minPart, maxPart, numParts;

};

#endif

