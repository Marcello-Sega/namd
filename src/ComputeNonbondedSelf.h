//-*-c++-*-
/***************************************************************************/
/*                                                                         */
/*              (C) Copyright 1996 The Board of Trustees of the            */
/*                          University of Illinois                         */
/*                           All Rights Reserved                           */
/*									   */
/***************************************************************************/

/***************************************************************************
 * DESCRIPTION: ComputeNonbondedSelf.h
 *
 ***************************************************************************/

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

