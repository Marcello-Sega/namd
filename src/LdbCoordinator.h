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

#include "ckdefs.h"
#include "chare.h"
#include "c++interface.h"

#include "converse.h"
#include "NamdTypes.h"

class PatchMap;
class ComputeMap;
class Sequencer;
class InitMsg;

enum {LDB_PATCHES = 1024};
enum {LDB_COMPUTES = 8192};

struct LdbStatsMsg : public comm_object
{
  int proc;
  int procLoad;
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

class LdbCoordinator : public groupmember
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
  void analyze(LdbStatsMsg *msg);
  void resume(LdbResumeMsg *msg);

private:
  int checkAndSendStats(void);
  void awakenSequencers(void);
  void printLocalLdbReport(void);
  void printLdbReport(void);

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
  int first_ldbcycle;

  LdbStatsMsg **statsMsgs;
};

#endif // LDBCOORDINATOR_H


/***************************************************************************
 * RCS INFORMATION:
 *
 *	$RCSfile $
 *	$Author $	$Locker:  $		$State: Exp $
 *	$Revision: 1.3 $	$Date: 1997/04/01 23:20:16 $
 *
 ***************************************************************************
 * REVISION HISTORY:
 *
 * $Log: LdbCoordinator.h,v $
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

