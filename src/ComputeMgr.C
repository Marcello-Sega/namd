/**
***  Copyright (c) 1995, 1996, 1997, 1998, 1999, 2000 by
***  The Board of Trustees of the University of Illinois.
***  All rights reserved.
**/

#include "charm++.h"

#include "ProcessorPrivate.h"

//#define DEBUGM
#define MIN_DEBUG_LEVEL 1
#include "Debug.h"

#include "BOCgroup.h"
#include "ComputeMgr.decl.h"
#include "ComputeMgr.h"
#include "ProxyMgr.decl.h"
#include "ProxyMgr.h"

#include "Node.h"
#include "ComputeMap.h"
#include "PatchMap.h"
#include "PatchMap.inl"

#include "Compute.h"
#include "ComputeNonbondedUtil.h"
#include "ComputeNonbondedSelf.h"
#include "ComputeNonbondedPair.h"
#include "ComputeAngles.h"
#include "ComputeDihedrals.h"
#include "ComputeImpropers.h"
#include "ComputeBonds.h"
#include "ComputeFullDirect.h"
#include "ComputeGlobal.h"
#include "ComputeGlobalMsgs.h"
#include "ComputeExt.h"
#include "ComputeDPMTA.h"
#include "ComputeDPME.h"
#include "ComputeDPMEMsgs.h"
#include "ComputePme.h"
#include "ComputeEField.h"
#include "ComputeSphericalBC.h"
#include "ComputeCylindricalBC.h"
#include "ComputeRestraints.h"
#include "ComputeConsForce.h"
#include "WorkDistrib.h"

/* include all of the specific masters we need here */
#include "FreeEnergyEnums.h"
#include "FreeEnergyAssert.h"
#include "FreeEnergyGroup.h"
#include "FreeEnergyVector.h"
#include "FreeEnergyRestrain.h"
#include "FreeEnergyRMgr.h"
#include "FreeEnergyLambda.h"
#include "FreeEnergyLambdMgr.h"

#include "GlobalMasterTest.h"
#include "GlobalMasterIMD.h"
#include "GlobalMasterTcl.h"
#include "GlobalMasterSMD.h"
#include "GlobalMasterEasy.h"
#include "GlobalMasterMisc.h"
#include "GlobalMasterFreeEnergy.h"

ComputeMgr::ComputeMgr()
{
  CpvAccess(BOCclass_group).computeMgr = thisgroup;
  computeGlobalObject = 0;
  computeDPMEObject = 0;
}

ComputeMgr::~ComputeMgr(void)
{
  ;
}

void ComputeMgr::updateComputes(int ep, CkGroupID chareID) {
  updateComputesReturnEP = ep;
  updateComputesReturnChareID = chareID;
  updateComputesCount = CkNumPes();

  if (CkMyPe()) { 
    iout << iPE << iERRORF << "updateComputes signaled on wrong Pe!\n" << endi;
    CkExit();
    return;
  }
  
#if CHARM_VERSION > 050402
  CkStartQD(CkIndex_ComputeMgr::updateComputes2((CkQdMsg*)0),&thishandle);
#else
  CkStartQD(CProxy_ComputeMgr::ckIdx_updateComputes2((CkQdMsg*)0),&thishandle);
#endif
}

void ComputeMgr::updateComputes2(CkQdMsg *msg) {
  delete msg;

  CProxy_WorkDistrib wd(CpvAccess(BOCclass_group).workDistrib);
  WorkDistrib  *workDistrib = wd.ckLocalBranch();
#if CHARM_VERSION > 050402
  workDistrib->saveComputeMapChanges(CkIndex_ComputeMgr::updateComputes3(),thisgroup);
#else
  workDistrib->saveComputeMapChanges(CProxy_ComputeMgr::ckIdx_updateComputes3(),thisgroup);
#endif
}

void ComputeMgr::updateComputes3() {
  CProxy_ComputeMgr(thisgroup).updateLocalComputes();
}

void ComputeMgr::updateLocalComputes() {
  ComputeMap *computeMap = ComputeMap::Object();
  CProxy_ProxyMgr pm(CpvAccess(BOCclass_group).proxyMgr);
  ProxyMgr *proxyMgr = pm.ckLocalBranch();

  computeFlag = new int[computeMap->numComputes()];

  for (int i=0; i<computeMap->numComputes(); i++) {
    DebugM(3, "updateLocalComputes("<<i<<") curnode="<<computeMap->node(i)
      <<" newnode="<<computeMap->newNode(i)<<"\n");
    computeFlag[i] = 0;
      
    if (computeMap->newNode(i) == CkMyPe() && computeMap->node(i) != CkMyPe()) {
      DebugM(4, "updateLocal - creating new computeID("<<i<<")\n");
      computeFlag[i] = 1;
      computeMap->setNode(i,computeMap->newNode(i));
      for (int n=0; n < computeMap->numPids(i); n++) {
	proxyMgr->createProxy(computeMap->pid(i,n));
      }
    } 
    else if (computeMap->node(i) == CkMyPe() && 
	(computeMap->newNode(i) != -1 && computeMap->newNode(i) != CkMyPe() )) {
      DebugM(4, "updateLocal - deleting computeID("<<i<<")\n");
      computeFlag[i] = -1;
      computeMap->setNode(i,computeMap->newNode(i));
    } else if (computeMap->newNode(i) != -1) {
      computeMap->setNode(i,computeMap->newNode(i));
    }
    computeMap->setNewNode(i,-1);
  }
 
  if (!CkMyPe()) {
#if CHARM_VERSION > 050402
      CkStartQD(CkIndex_ComputeMgr::updateLocalComputes2((CkQdMsg*)0), &thishandle);
#else
      CkStartQD(CProxy_ComputeMgr::ckIdx_updateLocalComputes2((CkQdMsg*)0), &thishandle);
#endif
  }
}

void
ComputeMgr::updateLocalComputes2(CkQdMsg *msg) {
  delete msg;
  CProxy_ComputeMgr(thisgroup).updateLocalComputes3();
}

void
ComputeMgr::updateLocalComputes3() {
  ComputeMap *computeMap = ComputeMap::Object();
  CProxy_ProxyMgr pm(CpvAccess(BOCclass_group).proxyMgr);
  ProxyMgr *proxyMgr = pm.ckLocalBranch();

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

  proxyMgr->removeUnusedProxies();

  DebugM(4, "msg to doneUpdateLocalComputes on Pe("<<CkMyPe()<<")\n");
  ComputeMap::Object()->checkMap();
  PatchMap::Object()->checkMap();

  if (!CkMyPe()) {
#if CHARM_VERSION > 050402
    CkStartQD(CkIndex_ComputeMgr::updateLocalComputes4((CkQdMsg*)0), &thishandle);
#else
    CkStartQD(CProxy_ComputeMgr::ckIdx_updateLocalComputes4((CkQdMsg*)0), &thishandle);
#endif
// added a new phase to build spanning tree after load balance
// was
//    CkStartQD(CProxy_ComputeMgr::ckIdx_doneUpdateLocalComputes(), &thishandle);
  }
  //CSendMsgBranch(ComputeMgr, doneUpdateLocalComputes, thisgroup, 0);
}

void
ComputeMgr::updateLocalComputes4(CkQdMsg *msg) {
  delete msg;
  CProxy_ComputeMgr(thisgroup).updateLocalComputes5();
}

void
ComputeMgr::updateLocalComputes5() {
  if (proxySendSpanning || proxyRecvSpanning )
    ProxyMgr::Object()->buildProxySpanningTree();
  if (!CkMyPe()) 
#if CHARM_VERSION > 050402
    CkStartQD(CkIndex_ComputeMgr::doneUpdateLocalComputes(), &thishandle);
#else
    CkStartQD(CProxy_ComputeMgr::ckIdx_doneUpdateLocalComputes(), &thishandle);
#endif
}

void ComputeMgr::doneUpdateLocalComputes() {

//  if (!--updateComputesCount) {
    DebugM(4, "doneUpdateLocalComputes on Pe("<<CkMyPe()<<")\n");
    void *msg = CkAllocMsg(0,0,0);
    CkSendMsgBranch(updateComputesReturnEP,msg,0,updateComputesReturnChareID);
//  }
}

//
void
ComputeMgr::createCompute(ComputeID i, ComputeMap *map)
{
    Compute *c;
    PatchID pid2[2];
    PatchIDList pids;
    int trans2[2];
    SimParameters *simParams = Node::Object()->simParameters;

    switch ( map->type(i) )
    {
      case computeNonbondedSelfType:
	c = new ComputeNonbondedSelf(i,map->computeData[i].pids[0].pid,
				     map->partition(i),map->partition(i)+1,
				     map->numPartitions(i)); // unknown delete
	++numNonbondedSelf;
	map->registerCompute(i,c);
	c->initialize();
	break;
      case computeNonbondedPairType:
	pid2[0] = map->computeData[i].pids[0].pid;
	trans2[0] = map->computeData[i].pids[0].trans;
	pid2[1] = map->computeData[i].pids[1].pid;
	trans2[1] = map->computeData[i].pids[1].trans;
	c = new ComputeNonbondedPair(i,pid2,trans2,
				     map->partition(i),map->partition(i)+1,
				     map->numPartitions(i)); // unknown delete
	++numNonbondedPair;
	map->registerCompute(i,c);
	c->initialize();
	break;
      case computeBondsType:
	PatchMap::Object()->homePatchIDList(pids);
	c = new ComputeBonds(i,pids); // unknown delete
	map->registerCompute(i,c);
	c->initialize();
	break;
      case computeAnglesType:
	PatchMap::Object()->homePatchIDList(pids);
	c = new ComputeAngles(i,pids); // unknown delete
	map->registerCompute(i,c);
	c->initialize();
	break;
      case computeDihedralsType:
	PatchMap::Object()->homePatchIDList(pids);
	c = new ComputeDihedrals(i,pids); // unknown delete
	map->registerCompute(i,c);
	c->initialize();
	break;
      case computeImpropersType:
	PatchMap::Object()->homePatchIDList(pids);
	c = new ComputeImpropers(i,pids); // unknown delete
	map->registerCompute(i,c);
	c->initialize();
	break;
      case computeSelfBondsType:
	c = new ComputeSelfBonds(i,map->computeData[i].pids[0].pid);
	map->registerCompute(i,c);
	c->initialize();
	break;
      case computeSelfAnglesType:
	c = new ComputeSelfAngles(i,map->computeData[i].pids[0].pid);
	map->registerCompute(i,c);
	c->initialize();
	break;
      case computeSelfDihedralsType:
	c = new ComputeSelfDihedrals(i,map->computeData[i].pids[0].pid);
	map->registerCompute(i,c);
	c->initialize();
	break;
      case computeSelfImpropersType:
	c = new ComputeSelfImpropers(i,map->computeData[i].pids[0].pid);
	map->registerCompute(i,c);
	c->initialize();
	break;
#ifdef DPMTA
      case computeDPMTAType:
	c = new ComputeDPMTA(i); // unknown delete
	map->registerCompute(i,c);
	c->initialize();
	break;
#endif
#ifdef DPME
      case computeDPMEType:
	c = computeDPMEObject = new ComputeDPME(i,this); // unknown delete
	map->registerCompute(i,c);
	c->initialize();
	break;
#endif
      case computePmeType:
	c = new ComputePme(i); // unknown delete
	map->registerCompute(i,c);
	c->initialize();
	break;
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
      case computeExtType:
	c = new ComputeExt(i); // unknown delete
	map->registerCompute(i,c);
	c->initialize();
	break;
      case computeEFieldType:
	c = new ComputeEField(i,map->computeData[i].pids[0].pid); // unknown delete
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
      case computeConsForceType:
        c = new ComputeConsForce(i,map->computeData[i].pids[0].pid);
        map->registerCompute(i,c);
        c->initialize();
        break;
      case computeConsTorqueType:
        c = new ComputeConsTorque(i,map->computeData[i].pids[0].pid);
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
  SimParameters *simParams = node->simParameters;
  int myNode = node->myid();

  numNonbondedSelf = 0;
  numNonbondedPair = 0;
  ComputeNonbondedUtil::select();

  if ( simParams->globalForcesOn && !myNode ) {
    DebugM(4,"Mgr running on Node "<<CkMyPe()<<"\n");
    /* create a master server to allow multiple masters */
    masterServerObject = new GlobalMasterServer(this,
    Node::Object()->workDistrib->getNumComputeGlobals());

    /* create the individual global masters */
    // masterServerObject->addClient(new GlobalMasterTest());
    if(simParams->tclForcesOn)
      masterServerObject->addClient(new GlobalMasterTcl());
    if(simParams->IMDon)
      masterServerObject->addClient(new GlobalMasterIMD());

    if(simParams->SMDOn)
      masterServerObject->addClient(
        new GlobalMasterSMD(simParams->SMDk, simParams->SMDVel,
			  simParams->SMDDir, simParams->SMDOutputFreq,
			  simParams->firstTimestep, simParams->SMDFile)
		  );
    if(simParams->miscForcesOn)
      masterServerObject->addClient(new GlobalMasterMisc());
    if ( simParams->freeEnergyOn )
      masterServerObject->addClient(new GlobalMasterFreeEnergy());
  }

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
      //      DebugM(1,"  pid " << map->computeData[i].pids[j] << '\n');
      if (!((j+1) % 6))
	DebugM(1,'\n');
    }
    DebugM(1,"\n---------------------------------------");
    DebugM(1,"---------------------------------------\n");

    createCompute(i, map);

  }

}


void ComputeMgr:: sendComputeGlobalConfig(ComputeGlobalConfigMsg *msg)
{
  (CProxy_ComputeMgr(CpvAccess(BOCclass_group).computeMgr)).recvComputeGlobalConfig(msg);
}

void ComputeMgr:: recvComputeGlobalConfig(ComputeGlobalConfigMsg *msg)
{
  if ( computeGlobalObject ) {
    computeGlobalObject->recvConfig(msg);
  }
  else if ( ! (PatchMap::Object())->numHomePatches() ) delete msg;
  else NAMD_die("ComputeMgr::computeGlobalObject is NULL!");
}

void ComputeMgr:: sendComputeGlobalData(ComputeGlobalDataMsg *msg)
{
  CProxy_ComputeMgr cm(CpvAccess(BOCclass_group).computeMgr);
#if CHARM_VERSION > 050402
  cm[0].recvComputeGlobalData(msg);
#else
  cm.recvComputeGlobalData(msg, 0);
#endif
}

void ComputeMgr:: recvComputeGlobalData(ComputeGlobalDataMsg *msg)
{
  if(masterServerObject) { // make sure it has been initialized
    masterServerObject->recvData(msg);
  } else NAMD_die("ComputeMgr::masterServerObject is NULL!");
}

void ComputeMgr:: sendComputeGlobalResults(ComputeGlobalResultsMsg *msg)
{
  (CProxy_ComputeMgr(CpvAccess(BOCclass_group).computeMgr)).recvComputeGlobalResults(msg);
}

void ComputeMgr:: recvComputeGlobalResults(ComputeGlobalResultsMsg *msg)
{
  if ( computeGlobalObject ) {
    computeGlobalObject->recvResults(msg);
  }
  else if ( ! (PatchMap::Object())->numHomePatches() ) delete msg;
  else NAMD_die("ComputeMgr::computeGlobalObject is NULL!");
}

void ComputeMgr:: sendComputeDPMEData(ComputeDPMEDataMsg *msg)
{
  if ( computeDPMEObject ) {
#ifdef DPME
    int node = computeDPMEObject->getMasterNode();
    CProxy_ComputeMgr cm(CpvAccess(BOCclass_group).computeMgr);
    cm.recvComputeDPMEData(msg,node);
#endif
  }
  else if ( ! (PatchMap::Object())->numHomePatches() ) delete msg;
  else NAMD_die("ComputeMgr::computeDPMEObject is NULL!");
}

void ComputeMgr:: recvComputeDPMEData(ComputeDPMEDataMsg *msg)
{
  if ( computeDPMEObject ) {
#ifdef DPME
    computeDPMEObject->recvData(msg);
#endif
  }
  else if ( ! (PatchMap::Object())->numHomePatches() ) delete msg;
  else NAMD_die("ComputeMgr::computeDPMEObject is NULL!");
}

void ComputeMgr:: sendComputeDPMEResults(ComputeDPMEResultsMsg *msg, int node)
{
  CProxy_ComputeMgr cm(CpvAccess(BOCclass_group).computeMgr);
#if CHARM_VERSION > 050402
  cm[node].recvComputeDPMEResults(msg);
#else
  cm.recvComputeDPMEResults(msg, node);
#endif
}

void ComputeMgr:: recvComputeDPMEResults(ComputeDPMEResultsMsg *msg)
{
  if ( computeDPMEObject ) {
#ifdef DPME
    computeDPMEObject->recvResults(msg);
#endif
  }
  else if ( ! (PatchMap::Object())->numHomePatches() ) delete msg;
  else NAMD_die("ComputeMgr::computeDPMEObject is NULL!");
}

#include "ComputeMgr.def.h"

