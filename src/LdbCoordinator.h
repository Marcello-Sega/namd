/**
***  Copyright (c) 1995, 1996, 1997, 1998, 1999, 2000 by
***  The Board of Trustees of the University of Illinois.
***  All rights reserved.
**/

#ifndef LDBCOORDINATOR_H
#define LDBCOORDINATOR_H

#include <stdio.h>

#include "charm++.h"

#include "converse.h"
#include "NamdTypes.h"
#include "BOCgroup.h"
#include "LdbCoordinator.decl.h"

class PatchMap;
class ComputeMap;
class Controller;
class Sequencer;
class computeInfo;
class patchInfo;
class processorInfo;

enum {LDB_PATCHES = 4096};
enum {LDB_COMPUTES = 16384};
enum {COMPUTEMAX = 16384};
enum {PATCHMAX = 4096};
enum {PROCESSORMAX = 512};


class LdbStatsMsg : public CMessage_LdbStatsMsg
{
public:
  int proc;
  double procLoad;
  int nPatches;
  int *pid;
  int *nAtoms;
  int nComputes;
  int *cid;
  float *computeTime;

  LdbStatsMsg(void);
  LdbStatsMsg(int npatches, int ncomputes);
  ~LdbStatsMsg(void);
  static void* pack(LdbStatsMsg* msg);
  static LdbStatsMsg* unpack(void *ptr);
};



struct LdbResumeMsg : public CMessage_LdbResumeMsg
{
  int dummy;
};

class LdbCoordinator : public BOCclass
{
public:
  LdbCoordinator();
  ~LdbCoordinator(void);
  static LdbCoordinator *Object()  { 
    return CpvAccess(LdbCoordinator_instance); 
  }

  void initialize(PatchMap *pmap, ComputeMap *cmap, int reinit=0);
  void patchLoad(PatchID id, int nAtoms, int timestep);
  void startWork(ComputeID id, int timestep);
  void endWork(ComputeID id, int timestep);
  void rebalance(Sequencer *seq, PatchID id);
  void rebalance(Controller *seq);
  void nodeDone(LdbResumeMsg *msg);
  void sendStats(LdbResumeMsg *msg);
  void analyze(LdbStatsMsg *msg);
  void updateComputesReady();
  void resume(LdbResumeMsg *msg);
  void resumeReady(CkQdMsg *msg);
  void resume2(LdbResumeMsg *msg);
  int steps(void) { return nLdbSteps; }

  // Public variables accessed by the idle-event functions
  double idleStart;
  double idleTime;

private:
  int checkAndGoToBarrier(void);
  void processStatistics(void);
  void awakenSequencers(void);
  int requiredProxies(PatchID id, int []);
  int buildData(void);
  void cleanUpData(void);
  void printRequiredProxies(PatchID id, FILE *fp);
  void printLocalLdbReport(void);
  void printLdbReport(const int nMoveableComputes);

  int stepsPerLdbCycle;
  int nLocalComputes;
  int nLocalPatches;
  int nPatchesReported;
  int nPatchesExpected;
  int nComputesReported;
  int nComputesExpected;
  int controllerReported;
  int controllerExpected;
  int nStatsMessagesReceived;
  int nStatsMessagesExpected;
  ComputeMap *computeMap;
  PatchMap *patchMap;
  int *patchNAtoms;
  Controller *controllerThread;
  Sequencer **sequencerThreads;
  double *computeStartTime;
  double *computeTotalTime;
  int ldbCycleNum;
  int nLdbSteps;
  int firstLdbStep;
  int nodesDone;

  LdbStatsMsg **statsMsgs;
  FILE *ldbStatsFP;
  double totalStartTime;
  double totalTime;
  computeInfo *computeArray;
  patchInfo *patchArray;
  processorInfo *processorArray;
};

#endif // LDBCOORDINATOR_H

