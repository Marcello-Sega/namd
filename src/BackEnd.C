/**
***  Copyright (c) 1995, 1996, 1997, 1998, 1999, 2000 by
***  The Board of Trustees of the University of Illinois.
***  All rights reserved.
**/

#include "BackEnd.h"
#include "ProcessorPrivate.h"
#include "common.h"
#include "Node.h"
#include "memusage.h"

#include <new>
#if defined(WIN32) && !defined(__CYGWIN__)
#include <new.h>
#endif

#ifdef USE_COMM_LIB
#include "ComlibManager.h"
#endif

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
#include "CollectionMaster.h"
#include "BroadcastMgr.decl.h"
#include "LdbCoordinator.decl.h"
#include "Sync.decl.h"

extern void _initCharm(int, char**);

float cpuTime_start;
float wallTime_start;

CkpvStaticDeclare(int,exitSchedHndlr);

extern "C" void exit_sched(void* msg)
{
  //  CmiPrintf("Exiting scheduler on %d\n",CmiMyPe());
  CsdExitScheduler();
}

static void register_exit_sched(void)
{
  CkpvInitialize(int,exitSchedHndlr);
  CkpvAccess(exitSchedHndlr) = CmiRegisterHandler((CmiHandler)exit_sched);
}

void BackEnd::ExitSchedOn(int pe)
{
  void* msg = CmiAlloc(CmiMsgHeaderSizeBytes);
  CmiSetHandler(msg,CkpvAccess(exitSchedHndlr));
  CmiSyncSendAndFree(pe,CmiMsgHeaderSizeBytes,(char *)msg);
}

#if defined(WIN32) && !defined(__CYGWIN__)
int NAMD_new_handler(size_t) {
#else
void NAMD_new_handler() {
#endif
  char tmp[100];
  sprintf(tmp,"Memory allocation failed on processor %d.",CmiMyPe());
  NAMD_die(tmp);
#if defined(WIN32) && !defined(__CYGWIN__)
  return 0;
#endif
}

// called on all procs
void all_init(int argc, char **argv)
{
#if defined(WIN32) && !defined(__CYGWIN__)
  _set_new_handler(NAMD_new_handler);
#else
  std::set_new_handler(NAMD_new_handler);
#endif
  ProcessorPrivateInit();
  register_exit_sched();
  _initCharm(argc, argv);  // message main Chare
}

// called on slave procs
void slave_init(int argc, char **argv)
{
  all_init(argc, argv);
  if (CkMyRank() < CkMyNodeSize()) 	// skip the communication thread
    CsdScheduler(-1);
}

// called by main on one or all procs
void BackEnd::init(int argc, char **argv) {
  ConverseInit(argc, argv, slave_init, 1, 1);  // calls slave_init on others
  cpuTime_start = CmiCpuTimer();
  wallTime_start = CmiWallTimer();
  if ( CmiMyPe() ) {
    slave_init(argc, argv);  // for procs that call main
    ConverseExit();  // should never return
  }
  all_init(argc, argv);

  // Create branch-office chares
  BOCgroup group;
  group.workDistrib = CProxy_WorkDistrib::ckNew();
  group.proxyMgr = CProxy_ProxyMgr::ckNew();
  group.patchMgr = CProxy_PatchMgr::ckNew();
  group.computeMgr = CProxy_ComputeMgr::ckNew();
  group.reductionMgr = CProxy_ReductionMgr::ckNew();
  group.computePmeMgr = CProxy_ComputePmeMgr::ckNew();
  group.computeExtMgr = CProxy_ComputeExtMgr::ckNew();
  group.sync = CProxy_Sync::ckNew();

  #if defined(NODEAWARE_PROXY_SPANNINGTREE) && defined(USE_NODEPATCHMGR)
  group.nodeProxyMgr = CProxy_NodeProxyMgr::ckNew();
  #endif 

#ifdef MEM_OPT_VERSION
  MasterHandlerInitMsg *initmsg8 = new MasterHandlerInitMsg;
  //initmsg8->master = collectionMaster;
  CkChareID collectionMasterHanlder = CProxy_CollectionMasterHandler::ckNew(initmsg8, 0);
#else
  CkChareID collectionMaster = CProxy_CollectionMaster::ckNew(0);
#endif

  SlaveInitMsg *initmsg7 = new SlaveInitMsg;
#ifndef MEM_OPT_VERSION
  initmsg7->master = collectionMaster;
#endif
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
  CmiPrintf("====================================================\n\n"
	    "WallClock: %f  CPUTime: %f  Memory: %f MB\n",
	    wallTime, cpuTime, memusage_MB());
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

