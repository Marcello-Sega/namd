/**
***  Copyright (c) 1995, 1996, 1997, 1998, 1999, 2000 by
***  The Board of Trustees of the University of Illinois.
***  All rights reserved.
**/

#include "charm++.h"

#include "ProcessorPrivate.h"

#define MIN_DEBUG_LEVEL 4
//#define DEBUGM
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
#include "ComputeNonbondedExcl.h"
#include "ComputeFullDirect.h"
#include "ComputeGlobal.h"
#include "ComputeGlobalMsgs.h"
#include "ComputeDPMTA.h"
#include "ComputeDPME.h"
#include "ComputeDPMEMsgs.h"
#include "ComputePme.h"
#include "ComputePmeMsgs.h"
#include "ComputeSphericalBC.h"
#include "ComputeCylindricalBC.h"
#include "ComputeRestraints.h"
#include "WorkDistrib.h"

ComputeMgr::ComputeMgr()
{
  CpvAccess(BOCclass_group).computeMgr = thisgroup;
  computeGlobalObject = 0;
  computeDPMEObject = 0;
  computePmeObject = 0;
}

ComputeMgr::~ComputeMgr(void)
{
  ;
}

void ComputeMgr::updateComputes(int ep, int chareID) {
  updateComputesReturnEP = ep;
  updateComputesReturnChareID = chareID;
  updateComputesCount = CkNumPes();

  if (CkMyPe()) { 
    iout << iPE << iERRORF << "updateComputes signaled on wrong Pe!\n" << endi;
    CkExit();
    return;
  }
  CkStartQD(CProxy_ComputeMgr::ckIdx_updateComputes2((CkQdMsg*)0),&thishandle);
}

void ComputeMgr::updateComputes2(CkQdMsg *msg) {
  delete msg;
  CProxy_WorkDistrib wd(CpvAccess(BOCclass_group).workDistrib);
  WorkDistrib  *workDistrib = wd.ckLocalBranch();
  workDistrib->saveComputeMapChanges(CProxy_ComputeMgr::ckIdx_updateComputes3(),thisgroup);
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
 
  DebugM(4, "updateComputes - totalComputes = "<<Compute::totalComputes<<"\n");
  if (!CkMyPe()) {
      CkStartQD(CProxy_ComputeMgr::ckIdx_updateLocalComputes2((CkQdMsg*)0), &thishandle);
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

  DebugM(4, "msg to doneUpdateLocalComputes on Pe("<<CkMyPe()<<")\n");
  ComputeMap::Object()->checkMap();
  PatchMap::Object()->checkMap();

  if (!CkMyPe()) {
    CkStartQD(CProxy_ComputeMgr::ckIdx_doneUpdateLocalComputes(), &thishandle);
  }
  //CSendMsgBranch(ComputeMgr, doneUpdateLocalComputes, thisgroup, 0);
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
      case computeNonbondedExclType:
	c = new ComputeNonbondedExcls(i); // unknown delete
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
	c = computePmeObject = new ComputePme(i,this); // unknown delete
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
  ComputeNonbondedUtil::select();

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


}


void ComputeMgr:: sendComputeGlobalConfig(ComputeGlobalConfigMsg *msg)
{
  CProxy_ComputeMgr(CpvAccess(BOCclass_group).computeMgr).recvComputeGlobalConfig(msg);
}

void ComputeMgr:: recvComputeGlobalConfig(ComputeGlobalConfigMsg *msg)
{
  if ( computeGlobalObject ) {
    computeGlobalObject->recvConfig(msg);
  }
  else if ( CkMyPe() >= (PatchMap::Object())->numPatches() ) delete msg;
  else NAMD_die("ComputeMgr::computeGlobalObject is NULL!");
}

void ComputeMgr:: sendComputeGlobalData(ComputeGlobalDataMsg *msg)
{
  CProxy_ComputeMgr cm(CpvAccess(BOCclass_group).computeMgr);
  cm.recvComputeGlobalData(msg, 0);
}

void ComputeMgr:: recvComputeGlobalData(ComputeGlobalDataMsg *msg)
{
  if ( computeGlobalObject ) {
    computeGlobalObject->recvData(msg);
  }
  else if ( CkMyPe() >= (PatchMap::Object())->numPatches() ) delete msg;
  else NAMD_die("ComputeMgr::computeGlobalObject is NULL!");
}

void ComputeMgr:: sendComputeGlobalResults(ComputeGlobalResultsMsg *msg)
{
  CProxy_ComputeMgr(CpvAccess(BOCclass_group).computeMgr).recvComputeGlobalResults(msg);
}

void ComputeMgr:: recvComputeGlobalResults(ComputeGlobalResultsMsg *msg)
{
  if ( computeGlobalObject ) {
    computeGlobalObject->recvResults(msg);
  }
  else if ( CkMyPe() >= (PatchMap::Object())->numPatches() ) delete msg;
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
  else if ( CkMyPe() >= (PatchMap::Object())->numPatches() ) delete msg;
  else NAMD_die("ComputeMgr::computeDPMEObject is NULL!");
}

void ComputeMgr:: recvComputeDPMEData(ComputeDPMEDataMsg *msg)
{
  if ( computeDPMEObject ) {
#ifdef DPME
    computeDPMEObject->recvData(msg);
#endif
  }
  else if ( CkMyPe() >= (PatchMap::Object())->numPatches() ) delete msg;
  else NAMD_die("ComputeMgr::computeDPMEObject is NULL!");
}

void ComputeMgr:: sendComputeDPMEResults(ComputeDPMEResultsMsg *msg, int node)
{
  CProxy_ComputeMgr cm(CpvAccess(BOCclass_group).computeMgr);
  cm.recvComputeDPMEResults(msg, node);
}

void ComputeMgr:: recvComputeDPMEResults(ComputeDPMEResultsMsg *msg)
{
  if ( computeDPMEObject ) {
#ifdef DPME
    computeDPMEObject->recvResults(msg);
#endif
  }
  else if ( CkMyPe() >= (PatchMap::Object())->numPatches() ) delete msg;
  else NAMD_die("ComputeMgr::computeDPMEObject is NULL!");
}

void ComputeMgr:: sendComputePmeData(ComputePmeDataMsg *msg)
{
  if ( computePmeObject ) {
    int node = computePmeObject->getMasterNode();
    CProxy_ComputeMgr cm(CpvAccess(BOCclass_group).computeMgr);
    cm.recvComputePmeData(msg,node);
  }
  else if ( CkMyPe() >= (PatchMap::Object())->numPatches() ) delete msg;
  else NAMD_die("ComputeMgr::computePmeObject is NULL!");
}

void ComputeMgr:: recvComputePmeData(ComputePmeDataMsg *msg)
{
  if ( computePmeObject ) {
    computePmeObject->recvData(msg);
  }
  else if ( CkMyPe() >= (PatchMap::Object())->numPatches() ) delete msg;
  else NAMD_die("ComputeMgr::computePmeObject is NULL!");
}

void ComputeMgr:: sendComputePmeResults(ComputePmeResultsMsg *msg, int node)
{
  CProxy_ComputeMgr cm(CpvAccess(BOCclass_group).computeMgr);
  cm.recvComputePmeResults(msg, node);
}

void ComputeMgr:: recvComputePmeResults(ComputePmeResultsMsg *msg)
{
  if ( computePmeObject ) {
    computePmeObject->recvResults(msg);
  }
  else if ( CkMyPe() >= (PatchMap::Object())->numPatches() ) delete msg;
  else NAMD_die("ComputeMgr::computePmeObject is NULL!");
}

#include "ComputeMgr.def.h"

