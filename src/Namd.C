/**
***  Copyright (c) 1995, 1996, 1997, 1998, 1999, 2000 by
***  The Board of Trustees of the University of Illinois.
***  All rights reserved.
**/

/*
   Namd() - launches BOC's Node, WorkDistrib, PatchMgr, ProxyMgr
   ReductionMgr, CollectionMgr
*/

#include "unistd.h"

#include "charm++.h"

#include "main.decl.h"
#include "main.h"
#include "BOCgroup.h"
#include "Namd.h"
#include "NamdState.h"
#include "Molecule.h"
#include "Parameters.h"
#include "SimParameters.h"
#include "ConfigList.h"
#include "Node.decl.h"
#include "Node.h"
#include "WorkDistrib.decl.h"
#include "WorkDistrib.h"
#include "PatchMgr.decl.h"
#include "PatchMgr.h"
#include "ComputeMgr.decl.h"
#include "ComputeMgr.h"
#include "ProxyMgr.decl.h"
#include "ProxyMgr.h"
#include "ReductionMgr.decl.h"
#include "ReductionMgr.h"
#include "CollectionMgr.decl.h"
#include "CollectionMgr.h"
#include "CollectionMaster.decl.h"
#include "CollectionMaster.h"
#include "BroadcastMgr.decl.h"
#include "BroadcastMgr.h"
#include "LdbCoordinator.decl.h"
#include "LdbCoordinator.h"

float Namd::cmiWallStart;
float Namd::cmiCpuStart;
int Namd::cmiFirstStart;
float Namd::cmiWallFirstStart;
float Namd::cmiCpuFirstStart;

// Namd(void ) is the constructor for the startup node.  It needs to
// read in file data,
Namd::Namd(void)
{
  namdState = new NamdState;

  BOCgroup group;

  // Create WorkDistrib and send it an empty message
  group.workDistrib = CProxy_WorkDistrib::ckNew();

  // Create ProxyMgr
  group.proxyMgr = CProxy_ProxyMgr::ckNew();

  // Create PatchMgr
  group.patchMgr = CProxy_PatchMgr::ckNew();

  // Create ComputeMgr
  group.computeMgr = CProxy_ComputeMgr::ckNew();

  // Create ReductionMgr
  group.reductionMgr = CProxy_ReductionMgr::ckNew();

  // Create Collection system
  CProxy_CollectionMaster coll(0);
  CkChareID collectionMaster = coll.ckGetChareId();
  SlaveInitMsg *initmsg7 = new SlaveInitMsg;
  initmsg7->master = collectionMaster;
  group.collectionMgr = CProxy_CollectionMgr::ckNew(initmsg7);

  // Create Broadcast system
  group.broadcastMgr = CProxy_BroadcastMgr::ckNew();

  // Create Load-balance coordinator
  group.ldbCoordinator = CProxy_LdbCoordinator::ckNew();

  // Create the Node object and send it the IDs of all the other
  // parallel objects.
  GroupInitMsg *msg = new GroupInitMsg;
  msg->group = group;

  // iout << iINFO << "Starting up nodes\n" << endi;
  nodeGroup = CProxy_Node::ckNew(msg);
}


// ~Namd(void) just needs to tell all the slave nodes to die.
Namd::~Namd(void)
{
  delete namdState;
  CkPrintf("Namd::~Namd() called\n");
}


// startup(char *) 
void Namd::startup(char *confFile)
{
  namdState->configFileInit(confFile);
  if (namdState->status()) {
    CkPrintf("Namd::startup() - could not initialize namdState from %s\n", 
      confFile);
    CkExit();
  }

  // Give our node[PE = 0] pointers to the data objects, so it can use them,
  // or send them on as messages elsewhere.
  Node::Object()->saveMolDataPointers(namdState);

  Node::messageStartUp(); // tell all nodes to startup
}

// last call of system
void Namd::namdDone(void) {
    Real CPUtime = CmiCpuTimer()-cmiCpuFirstStart;
    Real Walltime = CmiWallTimer()-cmiWallFirstStart;
    CkPrintf("==========================================\n"
	"WallClock : %f  CPUTime : %f \n",Walltime,CPUtime);
    CkExit();
}

