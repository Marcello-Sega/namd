/**
***  Copyright (c) 1995, 1996, 1997, 1998, 1999, 2000 by
***  The Board of Trustees of the University of Illinois.
***  All rights reserved.
**/

#ifndef COMPUTEPME_H
#define COMPUTEPME_H

#include "ComputeHomePatches.h"
#include "NamdTypes.h"

class ComputePmeDataMsg;
class ComputePmeResultsMsg;
class ComputePmeMaster;
class ComputeMgr;

class ComputePme : public ComputeHomePatches {
public:
  ComputePme(ComputeID c, ComputeMgr *m);
  virtual ~ComputePme();
  void doWork();
  void recvData(ComputePmeDataMsg *);
  void recvResults(ComputePmeResultsMsg *);

  ComputeMgr *comm;
  int getMasterNode(void) { return masterNode; }

 private:
  ComputePmeMaster *master;
  int masterNode;
  int numLocalAtoms;

};

#endif

