//-*-c++-*-
/***************************************************************************/
/*                                                                         */
/*              (C) Copyright 1996 The Board of Trustees of the            */
/*                          University of Illinois                         */
/*                           All Rights Reserved                           */
/*									   */
/***************************************************************************/

/***************************************************************************
 * DESCRIPTION:
 *
 ***************************************************************************/

#ifndef LDBCOORDINATOR_H
#define LDBCOORDINATOR_H

#include <stdio.h>

#include "ckdefs.h"
#include "chare.h"
#include "c++interface.h"

#include "converse.h"
#include "NamdTypes.h"
#include "elements.h"
#include "ComputeMgr.h"
#include "BOCgroup.h"
#include "Node.h"
#include "SimParameters.h"

class PatchMap;
class ComputeMap;
class Sequencer;
class InitMsg;

enum {LDB_PATCHES = 4096};
enum {LDB_COMPUTES = 16384};
enum {COMPUTEMAX = 16384};
enum {PATCHMAX = 4096};
enum {PROCESSORMAX = 512};


struct LdbStatsMsg : public comm_object
{
  int proc;
  double procLoad;
  int nPatches;
  int pid[LDB_PATCHES];
  int nAtoms[LDB_PATCHES];
  int nComputes;
  int cid[LDB_COMPUTES];
  float computeTime[LDB_COMPUTES];
};

struct LdbResumeMsg : public comm_object
{
  int dummy;
};

class LdbCoordinator : public BOCclass
{
public:
  LdbCoordinator(InitMsg *msg);
  ~LdbCoordinator(void);
  static LdbCoordinator *Object()  { return _instance; }

  void initialize(PatchMap *pmap, ComputeMap *cmap);
  void patchLoad(PatchID id, int nAtoms, int timestep);
  void startWork(ComputeID id, int timestep);
  void endWork(ComputeID id, int timestep);
  void rebalance(Sequencer *seq, PatchID id);
  void nodeDone(LdbResumeMsg *msg);
  void sendStats(LdbResumeMsg *msg);
  void analyze(LdbStatsMsg *msg);
  void updateComputesReady(DoneMsg *msg);
  void resume(LdbResumeMsg *msg);
  inline int balanceNow(int timestep);

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

  static LdbCoordinator *_instance;
  int stepsPerLdbCycle;
  int nLocalComputes;
  int nLocalPatches;
  int nPatchesReported;
  int nPatchesExpected;
  int nComputesReported;
  int nComputesExpected;
  int nStatsMessagesReceived;
  int nStatsMessagesExpected;
  ComputeMap *computeMap;
  PatchMap *patchMap;
  int *patchNAtoms;
  Sequencer **sequencerThreads;
  double *computeStartTime;
  double *computeTotalTime;
  Boolean ldbCycleNum;
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


inline int LdbCoordinator::balanceNow(int timestep)
{
  Node *node = Node::Object();
  const SimParameters *simParams = node->simParameters;
  const int numberOfSteps = simParams->N;
  int stepno = timestep - simParams->firstTimestep + 1;
  int firststep = firstLdbStep - simParams->firstTimestep;

  return 
    ( (node->numNodes() != 1) 
      && (simParams->ldbStrategy != LDBSTRAT_NONE) 
      && (stepno <= numberOfSteps)
      && (stepno >= firststep)
      && (
	  (stepno == firststep)
	  || (((stepno - 2*firststep) % simParams->ldbPeriod) == 0))
         );
}

#endif // LDBCOORDINATOR_H


/***************************************************************************
 * RCS INFORMATION:
 *
 *	$RCSfile $
 *	$Author $	$Locker:  $		$State: Exp $
 *	$Revision: 1.15 $	$Date: 1997/09/02 15:30:11 $
 *
 ***************************************************************************
 * REVISION HISTORY:
 *
 * $Log: LdbCoordinator.h,v $
 * Revision 1.15  1997/09/02 15:30:11  brunner
 * Changed some static arrays to allow only a max of 16384 computes
 *
 * Revision 1.14  1997/08/29 22:00:29  brunner
 * Load balancing improvements, and some diagnostic prints that need
 * to be removed, but won't occur if load balancing is off.
 *
 * Revision 1.13  1997/08/27 18:37:00  brunner
 * Load balancer end of computation sync added.  AlgSeven modified slightly.
 *
 * Revision 1.12  1997/07/08 15:48:10  milind
 * Made namd2 to work with Origin2000: Again...
 *
 * Revision 1.11  1997/04/16 23:44:03  brunner
 * Put ldbStrategy={none|refineonly|alg7}, ldbPeriod, and firstLdbStep
 * in SimParameters.
 *
 * Revision 1.10  1997/04/16 22:12:17  brunner
 * Fixed an LdbCoordinator bug, and cleaned up timing and Ldb output some.
 *
 * Revision 1.9  1997/04/10 09:14:00  ari
 * Final debugging for compute migration / proxy creation for load balancing.
 * Lots of debug code added, mostly turned off now.
 * Fixed bug in PositionBox when Patch had no dependencies.
 * Eliminated use of cout and misuse of iout in numerous places.
 *                                            Ari & Jim
 *
 * Revision 1.8  1997/04/08 23:00:01  brunner
 * Fixed problem with numComputes not equal to number of moveable computes
 *
 * Revision 1.7  1997/04/08 07:08:47  ari
 * Modification for dynamic loadbalancing - moving computes
 * Still bug in new computes or usage of proxies/homepatches.
 * Works if ldbStrategy is none as before.
 *
 * Revision 1.6  1997/04/07 00:04:07  brunner
 * LDB Alg7 added and computeMap modified on node 0
 *
 * Revision 1.5  1997/04/05 04:53:11  brunner
 * Rebalancer code linked in.  It is called, but the results are not currently
 * used.
 *
 * Revision 1.4  1997/04/04 17:31:42  brunner
 * New charm fixes for CommunicateConverse, and LdbCoordinator data file
 * output, required proxies, and idle time.
 *
 * Revision 1.3  1997/04/01 23:20:16  brunner
 * Collection on node 0 added
 *
 * Revision 1.2  1997/04/01 18:08:44  brunner
 * Made counts work right for first cycle
 *
 * Revision 1.1  1997/03/27 20:25:47  brunner
 * Changes for LdbCoordinator, the load balance control BOC
 *
 *
 ***************************************************************************/

