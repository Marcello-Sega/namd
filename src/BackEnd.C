/**
***  Copyright (c) 1995, 1996, 1997, 1998, 1999, 2000 by
***  The Board of Trustees of the University of Illinois.
***  All rights reserved.
**/

#include "BackEnd.h"

#include "converse.h"
#include "ProcessorPrivate.h"
#include "common.h"
#include "Node.h"

#include "charm++.h"
#include "main.decl.h"
#include "main.h"
#include "BOCgroup.h"
#include "WorkDistrib.decl.h"
#include "ProxyMgr.decl.h"
#include "PatchMgr.decl.h"
#include "ComputeMgr.decl.h"
#include "ReductionMgr.decl.h"
#include "CollectionMgr.decl.h"
#include "CollectionMaster.decl.h"
#include "CollectionMgr.h"
#include "BroadcastMgr.decl.h"
#include "LdbCoordinator.decl.h"

extern void _initCharm(int, char**);

float cpuTime_start;
float wallTime_start;

CpvStaticDeclare(int,exitSchedHndlr);

extern "C" void exit_sched(void* msg)
{
  //  CmiPrintf("Exiting scheduler on %d\n",CmiMyPe());
  CsdExitScheduler();
  CmiFree(msg);
}

static void register_exit_sched(void)
{
  CpvInitialize(int,exitSchedHndlr);
  CpvAccess(exitSchedHndlr) = CmiRegisterHandler((CmiHandler)exit_sched);
}

void BackEnd::ExitSchedOn(int pe)
{
  void* msg = CmiAlloc(CmiMsgHeaderSizeBytes);
  CmiSetHandler(msg,CpvAccess(exitSchedHndlr));
  CmiSyncSendAndFree(pe,CmiMsgHeaderSizeBytes,msg);
}

// called on all procs by namd_init()
void slave_init(int argc, char **argv)
{
  ProcessorPrivateInit();
  register_exit_sched();
  _initCharm(argc, argv);
  CsdScheduler(-1);
}

// called on all procs by front end
void BackEnd::init(int argc, char **argv) {
  ConverseInit(argc, argv, slave_init, 1, 1);
  cpuTime_start = CmiCpuTimer();
  wallTime_start = CmiWallTimer();
  if ( CmiMyPe() ) {
    slave_init(argc, argv);  // for procs that call main
    ConverseExit();  // should never return
  }
  ProcessorPrivateInit();
  register_exit_sched();
  _initCharm(argc, argv);  // message main Chare

  // Create branch-office chares
  BOCgroup group;
  group.workDistrib = CProxy_WorkDistrib::ckNew();
  group.proxyMgr = CProxy_ProxyMgr::ckNew();
  group.patchMgr = CProxy_PatchMgr::ckNew();
  group.computeMgr = CProxy_ComputeMgr::ckNew();
  group.reductionMgr = CProxy_ReductionMgr::ckNew();
  CProxy_CollectionMaster coll(0);
  CkChareID collectionMaster = coll.ckGetChareId();
  SlaveInitMsg *initmsg7 = new SlaveInitMsg;
  initmsg7->master = collectionMaster;
  group.collectionMgr = CProxy_CollectionMgr::ckNew(initmsg7);
  group.broadcastMgr = CProxy_BroadcastMgr::ckNew();
  group.ldbCoordinator = CProxy_LdbCoordinator::ckNew();
  GroupInitMsg *msg = new GroupInitMsg;
  msg->group = group;
  CProxy_Node::ckNew(msg);

}

// called on proc 0 by front end
void BackEnd::exit(void) {
  float cpuTime = CmiCpuTimer() - cpuTime_start;
  float wallTime = CmiWallTimer() - wallTime_start;
  CmiPrintf("==========================================\n"
  	"WallClock : %f  CPUTime : %f \n",wallTime,cpuTime);
  // char endmsg[1024];
  // sprintf(endmsg,"==========================================\n"
  // 	"WallClock : %f  CPUTime : %f \n",wallTime,cpuTime);
  // ConverseExit();  // Only kills this node.
  // CkExit(); CsdScheduler(-1);  // Kills everything?
  // CmiAbort("Exiting normally.\n");  // Actually works. 
  //                                     (Well, only on workstations -RKB)
  // 
  int i;
  for(i=1; i < CmiNumPes(); i++)
    ExitSchedOn(i);
  ConverseExit();
}

// start scheduler
void BackEnd::suspend(void) {
  CsdScheduler(-1);
}

// start quiescence detection to return to front end
void BackEnd::awaken(void) {
  Node::Object()->enableExitScheduler();
}

// start QD and scheduler
void BackEnd::barrier(void) {
  awaken();
  suspend();
}

