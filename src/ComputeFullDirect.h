/***************************************************************************/
/*      (C) Copyright 1996,1997 The Board of Trustees of the               */
/*                          University of Illinois                         */
/*                           All Rights Reserved                           */
/***************************************************************************/
/***************************************************************************
 * DESCRIPTION:  Computes electrostatics using efficient all atoms^2 calc
 *
 ***************************************************************************/

#ifndef COMPUTEFULLDIRECT_H
#define COMPUTEFULLDIRECT_H

#include "ComputeHomePatches.h"
#include "NamdTypes.h"

class SubmitReduction;

class ComputeFullDirect : public ComputeHomePatches {
private:
  SubmitReduction *reduction;
public:
  ComputeFullDirect(ComputeID c);
  virtual ~ComputeFullDirect();
  void doWork();
};

#endif

