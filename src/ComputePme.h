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
#include "MathArray.h"

class PmeRealSpace;
class ComputeMgr;
class SubmitReduction;
class PmeGridMsg;
#define PME_MAX_EVALS 3
typedef MathArray<double,7*PME_MAX_EVALS> PmeReduction;

class ComputePme : public ComputeHomePatches {
public:
  ComputePme(ComputeID c);
  virtual ~ComputePme();
  void doWork();
  void sendData(int, int*, int*, int*);
  void copyResults(PmeGridMsg *);
  void ungridForces();

 private:
  PmeGrid myGrid;
  int qsize, fsize, bsize;
  int fepOn, lesOn, lesFactor, numGrids;
  double **q_arr;
  char *f_arr;
  char *fz_arr;
  PmeReduction evir;
  SubmitReduction *reduction;
  int resultsRemaining;
  PmeRealSpace *myRealSpace[PME_MAX_EVALS];
  int numLocalAtoms;
  PmeParticle *localData;
  char *localPartition;
  int numGridAtoms[PME_MAX_EVALS];
  PmeParticle *localGridData[PME_MAX_EVALS];

};

#endif

