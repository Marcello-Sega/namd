/**
***  Copyright (c) 1995, 1996, 1997, 1998, 1999, 2000 by
***  The Board of Trustees of the University of Illinois.
***  All rights reserved.
**/

/*
   Forwards atoms to master node for force evaluation.
*/

#ifndef COMPUTEGLOBAL_H
#define COMPUTEGLOBAL_H

#include "ComputeHomePatches.h"
#include "NamdTypes.h"

class ComputeGlobalConfigMsg;
class ComputeGlobalDataMsg;
class ComputeGlobalResultsMsg;
class ComputeMgr;
class SubmitReduction;

class ComputeGlobal : public ComputeHomePatches {
public:
  ComputeGlobal(ComputeID, ComputeMgr*);
  virtual ~ComputeGlobal();
  void doWork();
  void recvConfig(ComputeGlobalConfigMsg *);
  void recvResults(ComputeGlobalResultsMsg *);

private:
  ComputeMgr *comm;

  void sendData();
  void configure(AtomIDList newaid, AtomIDList newgdef);

  AtomIDList aid;
  AtomIDList gdef;  // definitions of groups
  ResizeArray<BigReal> gmass;  // masses of groups

  int firsttime;
  SubmitReduction *reduction;
};

#endif

