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

#include "charm++.h"

#include "ProcessorPrivate.h"

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

#include "Compute.h"
#include "ComputeNonbondedUtil.h"
#include "ComputeNonbondedSelf.h"
#include "ComputeNonbondedPair.h"
#include "ComputeAngles.h"
#include "ComputeDihedrals.h"
#include "ComputeImpropers.h"
#include "ComputeBonds.h"
#include "ComputeNonbondedExcl.h"
#include "ComputeFullDirect.h"
#include "ComputeGlobal.h"
#include "ComputeGlobalMsgs.h"
#ifdef DPMTA
#include "ComputeDPMTA.h"
#endif
#include "ComputeSphericalBC.h"
#include "ComputeCylindricalBC.h"
#include "ComputeRestraints.h"
#include "ComputeSMD.h"
#include "WorkDistrib.h"

ComputeMgr::ComputeMgr(InitMsg *msg)
{
  delete msg;
  CpvAccess(BOCclass_group).computeMgr = thisgroup;
  computeGlobalObject = 0;
}

ComputeMgr::~ComputeMgr(void)
{
  ;
}

void ComputeMgr::updateComputes(int ep, int chareID) {
  updateComputesReturnEP = ep;
  updateComputesReturnChareID = chareID;
  updateComputesCount = CNumPes();

  if (CMyPe()) { 
    iout << iPE << iERRORF << "updateComputes signaled on wrong Pe!\n" << endi;
    CharmExit();
    return;
  }
  CStartQuiescence(GetEntryPtr(ComputeMgr,updateComputes2,QuiescenceMessage),thishandle);
}

void ComputeMgr::updateComputes2(QuiescenceMessage *msg) {
  delete msg;
  WorkDistrib  *workDistrib = CLocalBranch(WorkDistrib, CpvAccess(BOCclass_group).workDistrib);
  workDistrib->saveComputeMapChanges(
    GetEntryPtr(ComputeMgr, updateComputes3, DoneMsg), thisgroup
  );
}

void ComputeMgr::updateComputes3(DoneMsg *msg) {
  delete msg;
  RunMsg *runmsg = new (MsgIndex(RunMsg)) RunMsg;
  CBroadcastMsgBranch(ComputeMgr, updateLocalComputes, RunMsg, runmsg, thisgroup); 
}

void ComputeMgr::updateLocalComputes(RunMsg *msg) {
  delete msg;
  ComputeMap *computeMap = ComputeMap::Object();
  PatchMap *patchMap = PatchMap::Object();
  ProxyMgr *proxyMgr = CLocalBranch(ProxyMgr, CpvAccess(BOCclass_group).proxyMgr);

  computeFlag = new int[computeMap->numComputes()];

  for (int i=0; i<computeMap->numComputes(); i++) {
    DebugM(3, "updateLocalComputes("<<i<<") curnode="<<computeMap->node(i)
      <<" newnode="<<computeMap->newNode(i)<<"\n");
    computeFlag[i] = 0;
      
    if (computeMap->newNode(i) == CMyPe() && computeMap->node(i) != CMyPe()) {
      DebugM(4, "updateLocal - creating new computeID("<<i<<")\n");
      computeFlag[i] = 1;
      computeMap->setNode(i,computeMap->newNode(i));
      for (int n=0; n < computeMap->numPids(i); n++) {
	proxyMgr->createProxy(computeMap->pid(i,n));
      }
    } 
    else if (computeMap->node(i) == CMyPe() && 
	(computeMap->newNode(i) != -1 && computeMap->newNode(i) != CMyPe() )) {
      DebugM(4, "updateLocal - deleting computeID("<<i<<")\n");
      computeFlag[i] = -1;
      computeMap->setNode(i,computeMap->newNode(i));
    } else if (computeMap->newNode(i) != -1) {
      computeMap->setNode(i,computeMap->newNode(i));
    }
    computeMap->setNewNode(i,-1);
  }
 
  DebugM(4, "updateComputes - totalComputes = "<<Compute::totalComputes<<"\n");
  if (!CMyPe()) {
      CStartQuiescence(GetEntryPtr(ComputeMgr,updateLocalComputes2,QuiescenceMessage),thishandle);  
  }
}

void
ComputeMgr::updateLocalComputes2(QuiescenceMessage *msg) {
  delete msg;

  RunMsg *runmsg = new (MsgIndex(RunMsg)) RunMsg;
  CBroadcastMsgBranch(ComputeMgr, updateLocalComputes3, RunMsg, runmsg, thisgroup);
}

void
ComputeMgr::updateLocalComputes3(RunMsg *msg) {
  delete msg;

  ComputeMap *computeMap = ComputeMap::Object();

  for (int i=0; i<computeMap->numComputes(); i++) {
    if (1 == computeFlag[i]) {
      DebugM(4, "updateLocalCompute3() - create computeID(" << i << ")\n");
      createCompute(i, computeMap);
    }
    else if (-1 == computeFlag[i]) {
      // remove this compute
      DebugM(4, "updateLocalCompute3() - delete computeID(" << i << ")\n");
      delete computeMap->compute(i);
      computeMap->registerCompute(i,NULL);
    }
  }
  delete[] computeFlag;

  DebugM(4, "msg to doneUpdateLocalComputes on Pe("<<CMyPe()<<")\n");
  ComputeMap::Object()->checkMap();
  PatchMap::Object()->checkMap();

  if (!CMyPe()) {
    CStartQuiescence(GetEntryPtr(ComputeMgr,doneUpdateLocalComputes,DoneMsg),
		      thishandle);
  }
  //DoneMsg *donemsg = new (MsgIndex(DoneMsg)) DoneMsg;
  //CSendMsgBranch(ComputeMgr, doneUpdateLocalComputes, donemsg, thisgroup, 0);
}

void ComputeMgr::doneUpdateLocalComputes(DoneMsg *msg) {
  delete msg;
  

//  if (!--updateComputesCount) {
    DebugM(4, "doneUpdateLocalComputes on Pe("<<CMyPe()<<")\n");
    DoneMsg *donemsg = new (MsgIndex(DoneMsg)) DoneMsg;
    GeneralSendMsgBranch(updateComputesReturnEP, donemsg, 
      0, -1, updateComputesReturnChareID);
//  }
}

//
void
ComputeMgr::createCompute(ComputeID i, ComputeMap *map)
{
    Compute *c;
    PatchID pid2[2];
    int trans2[2];

    switch ( map->type(i) )
    {
      case computeNonbondedSelfType:
	c = new ComputeNonbondedSelf(i,map->computeData[i].pids[0].pid); // unknown delete
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
	++numNonbondedPair;
	map->registerCompute(i,c);
	c->initialize();
	break;
      case computeNonbondedExclType:
	c = new ComputeNonbondedExcls(i); // unknown delete
	map->registerCompute(i,c);
	c->initialize();
	break;
      case computeBondsType:
	c = new ComputeBonds(i); // unknown delete
	map->registerCompute(i,c);
	c->initialize();
	break;
      case computeAnglesType:
	c = new ComputeAngles(i); // unknown delete
	map->registerCompute(i,c);
	c->initialize();
	break;
      case computeDihedralsType:
	c = new ComputeDihedrals(i); // unknown delete
	map->registerCompute(i,c);
	c->initialize();
	break;
      case computeImpropersType:
	c = new ComputeImpropers(i); // unknown delete
	map->registerCompute(i,c);
	c->initialize();
	break;
#ifdef DPMTA
      case computeDPMTAType:
	c = new ComputeDPMTA(i); // unknown delete
	map->registerCompute(i,c);
	c->initialize();
	break;
      case computeDPMEType:
	// c = new ComputeDPME(i); // unknown delete
	// map->registerCompute(i,c);
	// c->initialize();
	break;
#endif
      case computeFullDirectType:
	c = new ComputeFullDirect(i); // unknown delete
	map->registerCompute(i,c);
	c->initialize();
	break;
      case computeGlobalType:
	c = computeGlobalObject = new ComputeGlobal(i,this); // unknown delete
	map->registerCompute(i,c);
	c->initialize();
	break;
      case computeSphericalBCType:
	c = new ComputeSphericalBC(i,map->computeData[i].pids[0].pid); // unknown delete
	map->registerCompute(i,c);
	c->initialize();
	break;
      case computeCylindricalBCType:
	c = new ComputeCylindricalBC(i,map->computeData[i].pids[0].pid); // unknown delete
	map->registerCompute(i,c);
	c->initialize();
	break;
      case computeRestraintsType:
	c = new ComputeRestraints(i,map->computeData[i].pids[0].pid); // unknown delete
	map->registerCompute(i,c);
	c->initialize();
	break;
      case computeSMDType:
	c = new ComputeSMD(i,map->computeData[i].pids[0].pid); // unknown delete
	map->registerCompute(i,c);
	c->initialize();
	break;
      
      default:
	break;
    }
}

   
void 
ComputeMgr::createComputes(ComputeMap *map)
{
  Node *node = Node::Object();
  int myNode = node->myid();

  numNonbondedSelf = 0;
  numNonbondedPair = 0;

  for(int i=0; i < map->nComputes; i++)
  {
    if ( ! ( i % 100 ) )
    {
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

  DebugM(4, "createComputes - total computes = "<<Compute::totalComputes<<"\n");

  // This doesn't really have to be here, but the output makes more sense.
  // It does have to happen after the molecule has been created.
  ComputeNonbondedUtil::select();

}


void ComputeMgr:: sendComputeGlobalConfig(ComputeGlobalConfigMsg *msg)
{
  CBroadcastMsgBranch(ComputeMgr, recvComputeGlobalConfig,ComputeGlobalConfigMsg, msg,
    CpvAccess(BOCclass_group).computeMgr);
}

void ComputeMgr:: recvComputeGlobalConfig(ComputeGlobalConfigMsg *msg)
{
  if ( computeGlobalObject ) {
    computeGlobalObject->recvConfig(msg);
  }
  else NAMD_die("ComputeMgr::computeGlobalObject is NULL!");
}

void ComputeMgr:: sendComputeGlobalData(ComputeGlobalDataMsg *msg)
{
  CSendMsgBranch(ComputeMgr, recvComputeGlobalData, ComputeGlobalDataMsg, msg,
    CpvAccess(BOCclass_group).computeMgr, 0);
}

void ComputeMgr:: recvComputeGlobalData(ComputeGlobalDataMsg *msg)
{
  if ( computeGlobalObject ) {
    computeGlobalObject->recvData(msg);
  }
  else NAMD_die("ComputeMgr::computeGlobalObject is NULL!");
}

void ComputeMgr:: sendComputeGlobalResults(ComputeGlobalResultsMsg *msg)
{
  CBroadcastMsgBranch(ComputeMgr, recvComputeGlobalResults, ComputeGlobalResultsMsg, msg,
    CpvAccess(BOCclass_group).computeMgr);
}

void ComputeMgr:: recvComputeGlobalResults(ComputeGlobalResultsMsg *msg)
{
  if ( computeGlobalObject ) {
    computeGlobalObject->recvResults(msg);
  }
  else NAMD_die("ComputeMgr::computeGlobalObject is NULL!");
}


#include "ComputeMgr.bot.h"


/***************************************************************************
 * RCS INFORMATION:
 *
 *	$RCSfile: ComputeMgr.C,v $
 *	$Author: brunner $	$Locker:  $		$State: Exp $
 *	$Revision: 1.1020 $	$Date: 1998/03/03 23:05:10 $
 *
 ***************************************************************************
 * REVISION HISTORY:
 *
 * $Log: ComputeMgr.C,v $
 * Revision 1.1020  1998/03/03 23:05:10  brunner
 * Changed include files for new simplified Charm++ include file structure.
 *
 * Revision 1.1019  1998/02/10 23:30:28  milind
 * Fixed to reflect the current changes to Charm++ translator.
 *
 * Revision 1.1018  1998/01/22 20:11:03  brunner
 * Modified the ComputeMap redistribution to send only new patch assignments.
 *
 * Revision 1.1017  1998/01/05 20:24:26  sergei
 * added #include "ComputeSMD.h";
 * added case computeSMDType in ComputeMgr::createCompute for SMD
 *
 * Revision 1.1016  1997/12/19 23:48:49  jim
 * Added Tcl interface for calculating forces.
 *
 * Revision 1.1015  1997/11/07 20:17:38  milind
 * Made NAMD to run on shared memory machines.
 *
 * Revision 1.1014  1997/10/01 16:46:48  milind
 * Removed old NAMD1 messaging and replaced it with new Message Streams library.
 *
 * Revision 1.1013  1997/04/22 04:25:58  jim
 * Added atomic restraints (harmonic constraints) via ComputeRestraints class.
 *
 * Revision 1.1012  1997/04/10 09:13:51  ari
 * Final debugging for compute migration / proxy creation for load balancing.
 * Lots of debug code added, mostly turned off now.
 * Fixed bug in PositionBox when Patch had no dependencies.
 * Eliminated use of cout and misuse of iout in numerous places.
 *                                            Ari & Jim
 *
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
