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

#include "ckdefs.h"
#include "chare.h"
#include "c++interface.h"

#define MIN_DEBUG_LEVEL 4
//#define DEBUGM
#include "Debug.h"

#include "BOCgroup.h"
#include "ComputeMgr.top.h"
#include "ComputeMgr.h"
#include "ProxyMgr.top.h"
#include "ProxyMgr.h"

#include "Node.h"
#include "ComputeMap.h"

#include "ComputeNonbondedUtil.h"
#include "ComputeNonbondedSelf.h"
#include "ComputeNonbondedPair.h"
#include "ComputeAngles.h"
#include "ComputeDihedrals.h"
#include "ComputeImpropers.h"
#include "ComputeBonds.h"
#include "ComputeNonbondedExcl.h"
#include "ComputeFullDirect.h"
#include "ComputeDPMTA.h"
#include "ComputeSphericalBC.h"
#include "ComputeCylindricalBC.h"
#include "WorkDistrib.h"

ComputeMgr::ComputeMgr(InitMsg *msg)
{
  delete msg;
  group.computeMgr = thisgroup;
}

ComputeMgr::~ComputeMgr(void)
{
  ;
}

void ComputeMgr::updateComputes(int ep, int chareID) {
  updateComputesReturnEP = ep;
  updateComputesReturnChareID = chareID;
  updateComputesCount = CNumPes();

  if (CMyPe()) { iout << iPE << iERRORF << 
    "updateComputes signaled on wrong Pe!\n" << endi;
    CharmExit();
    return;
  }
  WorkDistrib  *workDistrib = CLocalBranch(WorkDistrib, group.workDistrib);
  workDistrib->saveComputeMap(
    GetEntryPtr(ComputeMgr, updateComputes2), thisgroup
  );
}

void ComputeMgr::updateComputes2(DoneMsg *msg) {
  delete msg;
  CStartQuiescence(GetEntryPtr(ComputeMgr,updateComputes3),thishandle);
}

void ComputeMgr::updateComputes3(QuiescenceMessage *msg) {
  delete msg;
  DebugM(4, "Quiescence detected\n");
  RunMsg *runmsg = new (MsgIndex(RunMsg)) RunMsg;
  CBroadcastMsgBranch(ComputeMgr, updateLocalComputes, runmsg, thisgroup); 
  DebugM(4, "Broadcasting out to updateLocalComputes\n");
}

void ComputeMgr::updateLocalComputes(RunMsg *msg) {
  delete msg;
  ComputeMap *computeMap = ComputeMap::Object();
  PatchMap *patchMap = PatchMap::Object();
  ProxyMgr *proxyMgr = CLocalBranch(ProxyMgr, group.proxyMgr);

  computeFlag = new int[computeMap->numComputes()];

  DebugM(4, "updateLocalComputes() running\n");

  for (int i=0; i<computeMap->numComputes(); i++) {
    computeFlag[i] = 0;
    DebugM(4, " Compute#" << i << " node=" << computeMap->node(i) 
      << " newNode=" << computeMap->newNode(i) << "\n" );
      
    if (computeMap->newNode(i) == CMyPe() && computeMap->node(i) != CMyPe()) {
      DebugM(4, " Compute#" << i << " Add compute\n");
      computeFlag[i] = 1;
      computeMap->setNode(i,computeMap->newNode(i));
      for (int n=0; n < computeMap->numPids(i); n++) {
	proxyMgr->createProxy(computeMap->pid(i,n));
      }
    } 
    else if (computeMap->node(i) == CMyPe() && 
	(computeMap->newNode(i) != -1 && computeMap->newNode(i) != CMyPe() )) {
      DebugM(4, " Compute#" << i << " Delete compute\n");
      computeFlag[i] = -1;
      computeMap->setNode(i,computeMap->newNode(i));
    }
    computeMap->setNewNode(i,-1);
  }
 
  if (!CMyPe()) {
      CStartQuiescence(GetEntryPtr(ComputeMgr,updateLocalComputes2),thishandle);  }
  DebugM(4, "updateLocalComputes() first phase done\n");
}

void
ComputeMgr::updateLocalComputes2(QuiescenceMessage *msg) {
  delete msg;

  DebugM(4, "updateLocalComputes2() quiescence detected\n");
  RunMsg *runmsg = new (MsgIndex(RunMsg)) RunMsg;
  CBroadcastMsgBranch(ComputeMgr, updateLocalComputes3, runmsg, thisgroup);
}

void
ComputeMgr::updateLocalComputes3(RunMsg *msg) {
  delete msg;
  DebugM(4, "updateLocalComputes3() running\n");

  ComputeMap *computeMap = ComputeMap::Object();

  for (int i=0; i<computeMap->numComputes(); i++) {
    if (1 == computeFlag[i]) {
      DebugM(4, "trying to create compute #" << i << "\n");
      createCompute(i, computeMap);
    }
    else if (-1 == computeFlag[i]) {
      // remove this compute
      DebugM(4, "Removing compute #" << i 
	<< " Type=" << computeMap->type(i) << "\n");
      delete computeMap->compute(i);
      computeMap->registerCompute(i,NULL);
    }
  }
  delete[] computeFlag;

  DoneMsg *donemsg = new (MsgIndex(DoneMsg)) DoneMsg;
  CSendMsgBranch(ComputeMgr, doneUpdateLocalComputes, donemsg, thisgroup, 0);
}

void ComputeMgr::doneUpdateLocalComputes(DoneMsg *msg) {
  delete msg;

  DebugM(4, "doneUpdateLocalComputes() started\n");
  if (!--updateComputesCount) {
    DoneMsg *donemsg = new (MsgIndex(DoneMsg)) DoneMsg;
    GeneralSendMsgBranch(updateComputesReturnEP, donemsg, 
      0, -1, updateComputesReturnChareID);
  }
}


//
void
ComputeMgr::createCompute(ComputeID i, ComputeMap *map)
{
    Compute *c;
    PatchID pid2[2];
    int trans2[2];

  DebugM(2,"createComputes 2: looping " << i << "on type: " << map->computeData[i].type << "\n");
    switch ( map->type(i) )
    {
      case computeNonbondedSelfType:
	c = new ComputeNonbondedSelf(i,map->computeData[i].pids[0].pid); // unknown delete
	DebugM(4,"ComputeNonbondedSelf created.\n");
	++numNonbondedSelf;
	map->registerCompute(i,c);
	c->initialize();
	break;
      case computeNonbondedPairType:
	pid2[0] = map->computeData[i].pids[0].pid;
	trans2[0] = map->computeData[i].pids[0].trans;
	pid2[1] = map->computeData[i].pids[1].pid;
	trans2[1] = map->computeData[i].pids[1].trans;
	c = new ComputeNonbondedPair(i,pid2,trans2); // unknown delete
	DebugM(4,"ComputeNonbondedPair created.\n");
	++numNonbondedPair;
	map->registerCompute(i,c);
	c->initialize();
	break;
      case computeNonbondedExclType:
	c = new ComputeNonbondedExcls(i); // unknown delete
	DebugM(4,"ComputeNonbondedExcls created.\n");
	map->registerCompute(i,c);
	DebugM(3,"ComputeNonbondedExcls registered.\n");
	c->initialize();
	DebugM(3,"ComputeNonbondedExcls ready.\n");
	break;
      case computeBondsType:
	c = new ComputeBonds(i); // unknown delete
	DebugM(4,"ComputeBonds created.\n");
	map->registerCompute(i,c);
	c->initialize();
	break;
      case computeAnglesType:
	c = new ComputeAngles(i); // unknown delete
	DebugM(4,"ComputeAngles created.\n");
	map->registerCompute(i,c);
	c->initialize();
	break;
      case computeDihedralsType:
	c = new ComputeDihedrals(i); // unknown delete
	DebugM(4,"ComputeDihedrals created.\n");
	map->registerCompute(i,c);
	c->initialize();
	break;
      case computeImpropersType:
	c = new ComputeImpropers(i); // unknown delete
	DebugM(4,"ComputeImpropers created.\n");
	map->registerCompute(i,c);
	c->initialize();
	break;
#ifdef DPMTA
      case computeDPMTAType:
	c = new ComputeDPMTA(i); // unknown delete
	DebugM(4,"ComputeDPMTA created.\n");
	map->registerCompute(i,c);
	c->initialize();
	break;
      case computeDPMEType:
	// c = new ComputeDPME(i); // unknown delete
	DebugM(4,"ComputeDPME *NOT* created.\n");
	// map->registerCompute(i,c);
	// c->initialize();
	break;
#endif
      case computeFullDirectType:
	c = new ComputeFullDirect(i); // unknown delete
	DebugM(4,"ComputeFullDirect created.\n");
	map->registerCompute(i,c);
	c->initialize();
	break;
      case computeSphericalBCType:
	c = new ComputeSphericalBC(i,map->computeData[i].pids[0].pid); // unknown delete
	DebugM(4,"ComputeSphericalBC created.\n");
	map->registerCompute(i,c);
	c->initialize();
	break;
      case computeCylindricalBCType:
	c = new ComputeCylindricalBC(i,map->computeData[i].pids[0].pid); // unknown delete
	DebugM(4,"ComputeCylindricalBC created.\n");
	map->registerCompute(i,c);
	c->initialize();
	break;
      default:
	DebugM(10,"Unknown compute type not created!\n");
    }
}

   
void 
ComputeMgr::createComputes(ComputeMap *map)
{
  DebugM(2,"createComputes 0\n");
  Node *node = Node::Object();
  int myNode = node->myid();

  numNonbondedSelf = 0;
  numNonbondedPair = 0;

  DebugM(1,"---------------------------------------\n");
  DebugM(1,"---------------------------------------\n");

  DebugM(3,"nComputes = " << map->nComputes << '\n');
  DebugM(3,"nPatchBased = " << map->nPatchBased << '\n');
  DebugM(3,"nAtomBased = " << map->nAtomBased << '\n');
  DebugM(3,"nAllocated = " << map->nComputes << '\n');
  DebugM(2,"createComputes 1: looping " << map->nComputes << "\n");
  for(int i=0; i < map->nComputes; i++)
  {
    if ( ! ( i % 100 ) )
    {
      DebugM(4,"Created " << i << " compute objects so far.\n");
    }
    if ( map->computeData[i].node != myNode ) continue;
    DebugM(1,"Compute " << i << '\n');
    DebugM(1,"  node = " << map->computeData[i].node << '\n');
    DebugM(1,"  type = " << map->computeData[i].type << '\n');
    DebugM(1,"  patchBased = " << map->computeData[i].patchBased << '\n');
    DebugM(1,"  numPids = " << map->computeData[i].numPids << '\n');
    DebugM(1,"  numPidsAllocated = " << map->computeData[i].numPidsAllocated << '\n');
    for(int j=0; j < map->computeData[i].numPids; j++)
    {
      DebugM(1,"  pid " << map->computeData[i].pids[j] << '\n');
      if (!((j+1) % 6))
	DebugM(1,'\n');
    }
    DebugM(1,"\n---------------------------------------");
    DebugM(1,"---------------------------------------\n");

    createCompute(i, map);

  }
  DebugM(2,"createComputes 5: done looping\n");
  DebugM(4, numNonbondedSelf << " ComputeNonbondedSelf created\n");
  DebugM(4, numNonbondedPair << " ComputeNonbondedPair created\n");

  // This doesn't really have to be here, but the output makes more sense.
  // It does have to happen after the molecule has been created.
  ComputeNonbondedUtil::select();

}


#include "ComputeMgr.bot.h"


/***************************************************************************
 * RCS INFORMATION:
 *
 *	$RCSfile: ComputeMgr.C,v $
 *	$Author: ari $	$Locker:  $		$State: Exp $
 *	$Revision: 1.1011 $	$Date: 1997/04/08 07:08:16 $
 *
 ***************************************************************************
 * REVISION HISTORY:
 *
 * $Log: ComputeMgr.C,v $
 * Revision 1.1011  1997/04/08 07:08:16  ari
 * Modification for dynamic loadbalancing - moving computes
 * Still bug in new computes or usage of proxies/homepatches.
 * Works if ldbStrategy is none as before.
 *
 * Revision 1.1010  1997/03/27 20:25:41  brunner
 * Changes for LdbCoordinator, the load balance control BOC
 *
 * Revision 1.1009  1997/03/20 23:53:40  ari
 * Some changes for comments. Copyright date additions.
 * Hooks for base level update of Compute objects from ComputeMap
 * by ComputeMgr.  Useful for new compute migration functionality.
 *
 * Revision 1.1008  1997/03/19 11:54:09  ari
 * Add Broadcast mechanism.
 * Fixed RCS Log entries on files that did not have Log entries.
 * Added some register variables to Molecule and ComputeNonbondedExcl.C
 *
 ***************************************************************************/
