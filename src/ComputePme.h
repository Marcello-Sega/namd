/**
***  Copyright (c) 1995, 1996, 1997, 1998, 1999, 2000 by
***  The Board of Trustees of the University of Illinois.
***  All rights reserved.
**/

#ifndef COMPUTEPME_H
#define COMPUTEPME_H

#include "ComputeHomePatches.h"
#include "PmeBase.h"
#include "NamdTypes.h"

class PmeRealSpace;
class ComputeMgr;
class SubmitReduction;
class PmeGridMsg;
class PmeUntransMsg;

class ComputePme : public ComputeHomePatches {
public:
  ComputePme(ComputeID c);
  virtual ~ComputePme();
  void doWork();
  void sendData();
  void copyEnergy(PmeUntransMsg *);
  void copyResults(PmeGridMsg *);
  void ungridForces();

 private:
  PmeGrid myGrid;
  int qsize, fsize, bsize;
  double *q_arr;
  char *f_arr;
  double energy;
  double virial[6];
  SubmitReduction *reduction;
  int resultsRemaining;
  PmeRealSpace *myRealSpace;
  int numLocalAtoms;
  PmeParticle *localData;

};

#endif

