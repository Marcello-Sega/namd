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

class ComputePmeDataMsg;
class ComputePmeResultsMsg;
class ComputePmeMaster;
class PmeRealSpace;
class ComputeMgr;
class SubmitReduction;
class PmeUngridMsg;
class PmeUntransMsg;

class ComputePme : public ComputeHomePatches {
public:
  ComputePme(ComputeID c, ComputeMgr *m);
  virtual ~ComputePme();
  void doWork();
  void recvData(ComputePmeDataMsg *);
  void recvResults(ComputePmeResultsMsg *);
  void ComputePme::copyEnergy(PmeUntransMsg *);
  void ComputePme::copyResults(PmeUngridMsg *);
  void ComputePme::ungridForces();

  ComputeMgr *comm;
  int getMasterNode(void) { return masterNode; }

 private:
  ComputePmeMaster *master;
  int masterNode;

  PmeGrid myGrid;
  int qsize, bsize;
  double *q_arr;
  double energy;
  double virial[6];
  SubmitReduction *reduction;
  int resultsRemaining;
  PmeRealSpace *myRealSpace;
  int numLocalAtoms;
  PmeParticle *localData;

};

#endif

