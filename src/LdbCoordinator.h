/**
***  Copyright (c) 1995, 1996, 1997, 1998, 1999, 2000 by
***  The Board of Trustees of the University of Illinois.
***  All rights reserved.
**/

#ifndef LDBCOORDINATOR_H
#define LDBCOORDINATOR_H

#include <stdio.h>

#include <charm++.h>
#include <LBDatabase.h>

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

class LdbCoordinator : public BOCclass
{
public:
  LdbCoordinator();
  ~LdbCoordinator(void);
  static LdbCoordinator *Object()  { 
    return CpvAccess(LdbCoordinator_instance); 
  }

  void initialize(PatchMap *pmap, ComputeMap *cmap, int reinit=0);
  void createLoadBalancer();
  void patchLoad(PatchID id, int nAtoms, int timestep);
  void startWork(ComputeID id, int timestep);
  void endWork(ComputeID id, int timestep);
  void rebalance(Sequencer *seq, PatchID id);
  void rebalance(Controller *seq);
  void nodeDone(void);
  void updateComputesReady();
  void barrier(void);
  void resume(void);
  void resumeReady(CkQdMsg *msg);
  void resume2(void);
  int steps(void) { return nLdbSteps; }
  static void staticMigrateFn(LDObjHandle handle, int dest);
  static void staticStatsFn(LDOMHandle h, int state);
  static void staticQueryEstLoadFn(LDOMHandle h);
  static void staticReceiveAtSync(void* data);
  static void staticResumeFromSync(void* data);
  void ReceiveAtSync(void);
  void Migrate(LDObjHandle handle, int dest);
  void RecvMigrate(LdbMigrateMsg*);
  void ProcessMigrate(LdbMigrateMsg*);
  void ExpectMigrate(LdbMigrateMsg*);
  void ResumeFromSync(void);

private:
  struct Migration {
    int id;
    int from;
    int to;
    Migration* next;
  };

public:
  void ExecuteMigrations(void);
  void awakenSequencers(void);
  int requiredProxies(PatchID id, int []);
  void printRequiredProxies(PatchID id, FILE *fp);
  void printLocalLdbReport(void);

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

  int ldbCycleNum;
  int nLdbSteps;
  int firstLdbStep;
  int takingLdbData;
  int nodesDone;

  FILE *ldbStatsFP;
  computeInfo *computeArray;
  patchInfo *patchArray;
  processorInfo *processorArray;
  LBDatabase *theLbdb;
  LDOMid myOMid;
  LDOMHandle myHandle;
  LDObjHandle* objHandles;
  int nRegisteredObjs;
  LDBarrierClient ldBarrierHandle;
  int reg_all_objs;
  LDObjHandle* patchHandles;
  Migration* migrations;
};

class LdbMigrateMsg : public CMessage_LdbMigrateMsg
{
public:
  LDObjHandle handle;
  int from;
  int to;
};


#endif // LDBCOORDINATOR_H

