/**
***  Copyright (c) 1995, 1996, 1997, 1998, 1999, 2000 by
***  The Board of Trustees of the University of Illinois.
***  All rights reserved.
**/

#include <stdlib.h>
#include <charm++.h>
#include "NamdCentLB.h"
#include "NamdNborLB.h"

#include "InfoStream.h"

#include "HomePatch.h"
#include "LdbCoordinator.decl.h"
#include "LdbCoordinator.h"
#include "NamdTypes.h"
#include "Node.h"
#include "SimParameters.h"
#include "PatchMap.inl"
#include "ComputeMap.h"
//#define DEBUGM
#define MIN_DEBUG_LEVEL 3
#include "Debug.h"
#include "Controller.h"
#include "Sequencer.h"
#include "RefineOnly.h"
#include "ComputeMgr.h"
#include "packmsg.h"

#include "elements.h"
#include "ComputeMgr.decl.h"

#define DEBUG_LEVEL 4

 
void LdbCoordinator::staticMigrateFn(LDObjHandle handle, int dest)
{
   ((LdbCoordinator*)handle.omhandle.user_ptr)->Migrate(handle,dest);
}

void LdbCoordinator::Migrate(LDObjHandle handle, int dest)
{
  LdbMigrateMsg* msg = new LdbMigrateMsg;
  msg->handle = handle;
  msg->from = CkMyPe();
  msg->to = dest;
#if CHARM_VERSION > 050402
  CProxy_LdbCoordinator ldbProxy(thisgroup);
  ldbProxy[CkMyPe()].RecvMigrate(msg);
#else
  CProxy_LdbCoordinator(thisgroup).RecvMigrate(msg,CkMyPe());
#endif
}

void LdbCoordinator::staticStatsFn(LDOMHandle h, int state)
{
  CkPrintf("I'm supposed to set stats\n");
}

void LdbCoordinator::staticQueryEstLoadFn(LDOMHandle h)
{
  CkPrintf("I'm supposed to query load\n");
}

void LdbCoordinator::staticReceiveAtSync(void* data)
{
  ((LdbCoordinator*)data)->ReceiveAtSync();
}

void LdbCoordinator::ReceiveAtSync()
{
  theLbdb->RegisteringObjects(myHandle);
}

void LdbCoordinator::staticResumeFromSync(void* data)
{
  ((LdbCoordinator*)data)->ResumeFromSync();
}

void LdbCoordinator::ResumeFromSync()
{
  theLbdb->DoneRegisteringObjects(myHandle);
#if CMK_PERSISTENT_COMM
  if (takingLdbData) {
//CmiPrintf("[%d] CmiDestoryPersistent\n", CkMyPe());

    HomePatchList *hpl = PatchMap::Object()->homePatchList();
    ResizeArrayIter<HomePatchElem> ai(*hpl);
    for (ai=ai.begin(); ai != ai.end(); ai++) {
      HomePatch *patch = (*ai).patch;
      patch->destoryPersistComm();
    }
  }
#endif
  CProxy_LdbCoordinator cl(thisgroup);
#if CHARM_VERSION > 050402
  cl[0].nodeDone();
#else
  cl.nodeDone(0);
#endif
}

LdbCoordinator::LdbCoordinator()
{
  if (CpvAccess(LdbCoordinator_instance) == NULL) {
    CpvAccess(LdbCoordinator_instance) = this;
  } else {
    iout << iFILE << iERROR << iPE 
	 << "LdbCoordinator instanced twice on same node!" << endi;
    CkExit();
  }
  
#if 0
  // Create a load balancer
  if (CkMyPe() == 0) {
    //   CreateCentralLB();
    CreateNamdCentLB();
    //   CreateNamdNborLB();
  }
#endif

  ldbCycleNum = 1;
  takingLdbData = 1;
  nLocalComputes = nLocalPatches = 0;
  patchNAtoms = (int *) NULL;
  sequencerThreads = (Sequencer **) NULL;
  ldbStatsFP = NULL;
  computeArray = NULL;
  patchArray = NULL;
  processorArray = NULL;

  nodesDone = 0;

  // Register self as an object manager for new charm++ balancer framework
  theLbdb = CProxy_LBDatabase::ckLocalBranch(lbdb);
#if CHARM_VERSION > 050403
  myOMid.id.idx = 1;
#else
  myOMid.id = 1;
#endif
  LDCallbacks cb = { (LDMigrateFn)staticMigrateFn,
		     (LDStatsFn)staticStatsFn,
		     (LDQueryEstLoadFn)staticQueryEstLoadFn
                   };
  myHandle = theLbdb->RegisterOM(myOMid,(void*)this,cb);

  // Add myself as a local barrier receiver, so I know when I might
  // be registering objects.
  theLbdb->AddLocalBarrierReceiver((LDBarrierFn)staticReceiveAtSync,
				   (void*)this);;

  // Also, add a local barrier client, to trigger load balancing
  ldBarrierHandle = theLbdb->
    AddLocalBarrierClient((LDResumeFn)staticResumeFromSync,
			  (void*)this);
  objHandles = 0;
  reg_all_objs = 1;
  migrations = 0;
}

LdbCoordinator::~LdbCoordinator(void)
{
  delete [] patchNAtoms;
  delete [] sequencerThreads;
  delete [] objHandles;
  if (CkMyPe() == 0)
  {
    delete [] computeArray;
    delete [] patchArray;
    delete [] processorArray;
  }
  if (ldbStatsFP)
    fclose(ldbStatsFP);

}

void LdbCoordinator::createLoadBalancer()
{
  const SimParameters *simParams = Node::Object()->simParameters;
  if (simParams->ldbStrategy == LDBSTRAT_ALGNBOR) 
    CreateNamdNborLB();
  else {
    //   CreateCentralLB();
    CreateNamdCentLB();
  }
}

void LdbCoordinator::initialize(PatchMap *pMap, ComputeMap *cMap, int reinit)
{
  const SimParameters *simParams = Node::Object()->simParameters;

#if 0
  static int lbcreated = 0;
  // PE0 first time Create a load balancer
  if (CkMyPe() == 0 && !lbcreated) {
    if (simParams->ldbStrategy == LDBSTRAT_ALGNBOR) 
      CreateNamdNborLB();
    else {
      //   CreateCentralLB();
      CreateNamdCentLB();
    }
    lbcreated = 1;
  }
#endif

  //  DebugM(10,"stepsPerLdbCycle initialized\n");
  stepsPerLdbCycle = simParams->ldbPeriod;
  firstLdbStep = simParams->firstLdbStep;

  computeMap = cMap;
  patchMap = pMap;

  // Set the number of received messages correctly for node 0

  nStatsMessagesExpected = Node::Object()->numNodes();
  nStatsMessagesReceived = 0;

  delete [] patchNAtoms;  // Depends on delete NULL to do nothing
  patchNAtoms = new int[pMap->numPatches()];

  typedef Sequencer *seqPtr;

  if ( ! reinit ) {
    delete [] sequencerThreads;  // Depends on delete NULL to do nothing
    sequencerThreads = new seqPtr[pMap->numPatches()];
  }

  nLocalPatches=0;

  int i;
  for(i=0;i<pMap->numPatches();i++)
  {
    if (pMap->node(i) == Node::Object()->myid())
    {
      nLocalPatches++;
      patchNAtoms[i]=0;
    } else {
      patchNAtoms[i]=-1;
    }
    if ( ! reinit ) sequencerThreads[i]=NULL;
  }
  if ( ! reinit ) controllerThread = NULL;
  if (nLocalPatches != pMap->numHomePatches())
    NAMD_die("Disaggreement in patchMap data.\n");
 
  nLocalComputes = 0;
  for(i=0;i<cMap->numComputes();i++)  {
    if ( (cMap->node(i) == Node::Object()->myid())
	 && ( (cMap->type(i) == computeNonbondedPairType)
	      || (cMap->type(i) == computeSelfBondsType)
	      || (cMap->type(i) == computeSelfAnglesType)
	      || (cMap->type(i) == computeSelfDihedralsType)
	      || (cMap->type(i) == computeSelfImpropersType)
	      || (cMap->type(i) == computeNonbondedSelfType) ) ) {
      nLocalComputes++;
    }
  }
  
  // New LB frameworks registration

  // Allocate data structure to save incoming migrations.  Processor
  // zero will get all migrations

  // If this is the first time through, we need it register patches
  if (reg_all_objs) {
    // Tell the lbdb that I'm registering objects, until I'm done
    // registering them.
    theLbdb->RegisteringObjects(myHandle);
    
    patchHandles = new LDObjHandle[nLocalPatches];
    int patch_count=0;
    int i;
    for(i=0;i<pMap->numPatches();i++)
      if (pMap->node(i) == Node::Object()->myid()) {
	LDObjid elemID;
	elemID.id[0] = i;
	elemID.id[1] = elemID.id[2] = elemID.id[3] = -2;

	if (patch_count >= nLocalPatches) {
	  iout << iFILE << iERROR << iPE 
	       << "LdbCoordinator found too many local patches!" << endi;
	  CkExit();
	}
	patchHandles[patch_count] 
	  = theLbdb->RegisterObj(myHandle,elemID,0,0);
	patch_count++;
      }
  
    // Allocate new object handles
    if (objHandles == 0) {
      objHandles = new LDObjHandle[cMap->numComputes()];
      for(i=0;i<cMap->numComputes();i++)
	objHandles[i].id.id[0] = -1; // Use -1 to mark unused entries

      // Register computes
      for(i=0;i<cMap->numComputes();i++)  {
	if ( (cMap->node(i) == Node::Object()->myid())
	     && ( (cMap->type(i) == computeNonbondedPairType)
	          || (cMap->type(i) == computeSelfBondsType)
	          || (cMap->type(i) == computeSelfAnglesType)
	          || (cMap->type(i) == computeSelfDihedralsType)
	          || (cMap->type(i) == computeSelfImpropersType)
		  || (cMap->type(i) == computeNonbondedSelfType) ) ) {
	  // Register the object with the load balancer
	  // Store the depended patch IDs in the rest of the element ID
	  LDObjid elemID;
	  elemID.id[0] = i;
	
	  if (cMap->numPids(i) > 2)
	    elemID.id[3] = cMap->pid(i,2);
	  else elemID.id[3] = -1;

	  if (cMap->numPids(i) > 1)
	    elemID.id[2] =  cMap->pid(i,1);
	  else elemID.id[2] = -1;

	  if (cMap->numPids(i) > 0)
	    elemID.id[1] =  cMap->pid(i,0);
	  else elemID.id[1] = -1;

	  objHandles[i] = theLbdb->RegisterObj(myHandle,elemID,0,1);
	}
      }
    }
    theLbdb->DoneRegisteringObjects(myHandle);
    reg_all_objs = 0;
  }

  // Fixup to take care of the extra timestep at startup
  // This is pretty ugly here, but it makes the count correct
  if (ldbCycleNum==1)
  {
    nLdbSteps = firstLdbStep;
    takingLdbData = 0;
    theLbdb->CollectStatsOff();
  }
  else if ((ldbCycleNum<=4) || (! takingLdbData))
  {
    nLdbSteps = firstLdbStep;
    takingLdbData = 1;
    theLbdb->CollectStatsOn();
  }
  else 
  {
    nLdbSteps = stepsPerLdbCycle - firstLdbStep;
    takingLdbData = 0;
    theLbdb->CollectStatsOff();
  }

  // in the case when trace is off at the beginning,
  // only turn trace of from after the first LB to the firstLdbStep after 
  // the second LB.
  // 1    2   3     4     5          6   7
  // off on Alg7 refine refine ...  on
#if CHARM_VERSION >= 050606
  if (traceAvailable()) {
    static int specialTracing = 0;
    if (ldbCycleNum == 1 && traceIsOn() == 0)  specialTracing = 1;
    if (specialTracing) {
      if (ldbCycleNum == 3) traceBegin();
      if (ldbCycleNum == 6) traceEnd();
    }
  }
#endif

  nPatchesReported = 0;
  nPatchesExpected = nLocalPatches;
  nComputesReported = 0;
  nComputesExpected = nLocalComputes * nLdbSteps;
  controllerReported = 0;
  controllerExpected = ! CkMyPe();

  if (CkMyPe() == 0)
  {
    if (computeArray == NULL)
      computeArray = new computeInfo[computeMap->numComputes()];
    if (patchArray == NULL)
      patchArray = new patchInfo[patchMap->numPatches()];
    if (processorArray == NULL)
      processorArray = new processorInfo[CkNumPes()];
  }
    
  theLbdb->ClearLoads();
}

void LdbCoordinator::patchLoad(PatchID id, int nAtoms, int /* timestep */)
{
  if (patchNAtoms[id] != -1) {
    patchNAtoms[id] = nAtoms;
    nPatchesReported++;
  } else {
    DebugM(10, "::patchLoad() Unexpected patch reporting in\n");
  }
}

void LdbCoordinator::startWork(ComputeID id, int /* timestep */ )
{
  theLbdb->ObjectStart(objHandles[id]);
}

void LdbCoordinator::endWork(ComputeID id, int /* timestep */)
{
  theLbdb->ObjectStop(objHandles[id]);
  nComputesReported++;
}

void LdbCoordinator::rebalance(Sequencer *seq, PatchID pid)
{
  if (Node::Object()->simParameters->ldbStrategy == LDBSTRAT_NONE)
    return;

  sequencerThreads[pid] = seq;
  seq->suspend();
}

void LdbCoordinator::rebalance(Controller *c)
{
  if (Node::Object()->simParameters->ldbStrategy == LDBSTRAT_NONE)
    return;

  DebugM(3, "Controller reached load balance barrier.\n");
  controllerReported = 1;
  controllerThread = c;

  CProxy_LdbCoordinator(thisgroup).barrier();

  CthSuspend();
}

void LdbCoordinator::barrier(void)
{
  if ( (nPatchesReported != nPatchesExpected) 
       || (nComputesReported != nComputesExpected)
       || (controllerReported != controllerExpected) )
  {
    NAMD_bug("Load balancer received wrong number of events.\n");
  }

  theLbdb->AtLocalBarrier(ldBarrierHandle);
}

void LdbCoordinator::nodeDone(void)
{
  nodesDone++;

  if (nodesDone==Node::Object()->numNodes()) {
    nodesDone=0;
    ExecuteMigrations();
  }
}

void LdbCoordinator::ExecuteMigrations(void)
{

  while (migrations) {
    computeMap->setNewNode(migrations->id,migrations->to);
    Migration *const next = migrations->next;
    delete migrations;
    migrations = next;
  }
 
 // computeMgr->updateComputes() call only on Node(0) i.e. right here
  // This will barrier for all Nodes - (i.e. Computes must be
  // here and with proxies before anyone can start up

  CProxy_ComputeMgr cm(CpvAccess(BOCclass_group).computeMgr);
  ComputeMgr *computeMgr = cm.ckLocalBranch();
#if CHARM_VERSION > 050402
  computeMgr->updateComputes(CkIndex_LdbCoordinator::
                             updateComputesReady(),thisgroup);
#else
  computeMgr->updateComputes(CProxy_LdbCoordinator::
			     ckIdx_updateComputesReady(),thisgroup);
#endif
}

void LdbCoordinator::RecvMigrate(LdbMigrateMsg* m)
{
  // This method receives the migration from the framework,
  // unregisters it, and
  // sends it to PE 0, which is where NAMD coordinates all of the
  // load balancing
  const int id = m->handle.id.id[0];

  theLbdb->UnregisterObj(objHandles[id]);
  objHandles[id].id.id[0] = -1;

#if CHARM_VERSION > 050402
  CProxy_LdbCoordinator ldbProxy(thisgroup);
  ldbProxy[0].ProcessMigrate(m);
#else
  CProxy_LdbCoordinator(thisgroup).ProcessMigrate(m,0);
#endif
}

void LdbCoordinator::ProcessMigrate(LdbMigrateMsg* m)
{
  // On PE 0, we store the migration, and then forward it on to the
  // destination PE

  Migration* new_migration = new Migration;

  new_migration->id = m->handle.id.id[0];
  new_migration->from = m->from;
  new_migration->to = m->to;
  new_migration->next = migrations;
  migrations = new_migration;

#if CHARM_VERSION > 050402
  CProxy_LdbCoordinator  ldbProxy(thisgroup);
  ldbProxy[m->to].ExpectMigrate(m);
#else
  CProxy_LdbCoordinator(thisgroup).ExpectMigrate(m,m->to);
#endif
}

void LdbCoordinator::ExpectMigrate(LdbMigrateMsg* m)
{
  objHandles[m->handle.id.id[0]] 
    = theLbdb->RegisterObj(myHandle,m->handle.id,0,1);

  theLbdb->Migrated(objHandles[m->handle.id.id[0]]);

  delete m;
}

void LdbCoordinator::updateComputesReady() {
  DebugM(3,"updateComputesReady()\n");

  CProxy_LdbCoordinator(thisgroup).resume();
#if CHARM_VERSION > 050402
  CkStartQD(CkIndex_LdbCoordinator::resumeReady((CkQdMsg*)0),&thishandle);
#else
  CkStartQD(CProxy_LdbCoordinator::ckIdx_resumeReady((CkQdMsg*)0),&thishandle);
#endif
}

void LdbCoordinator::resume(void)
{
  DebugM(3,"resume()\n");

  //  printLocalLdbReport();

  ldbCycleNum++;
  initialize(PatchMap::Object(),ComputeMap::Object(),1);
}

void LdbCoordinator::resumeReady(CkQdMsg *msg) {

  DebugM(3,"resumeReady()\n");
  delete msg;

  CProxy_LdbCoordinator(thisgroup).resume2();
}

void LdbCoordinator::resume2(void)
{
  DebugM(3,"resume2()\n");

  awakenSequencers();
}

void LdbCoordinator::awakenSequencers()
{
  if (controllerThread)
  {
    controllerThread->awaken();
    controllerThread = NULL;
  }
  for(int i=0; i < patchMap->numPatches(); i++)
  {
    if (sequencerThreads[i])
    {
      sequencerThreads[i]->awaken();
    }
    sequencerThreads[i]= NULL;
  }
}

// Figure out which proxies we will definitely create on other
// nodes, without regard for non-bonded computes.  This code is swiped
// from ProxyMgr, and changes there probable need to be propagated here.

int LdbCoordinator::requiredProxies(PatchID id, int neighborNodes[])
{
  enum proxyHere { No, Yes };
  int numNodes = CkNumPes();
  proxyHere *proxyNodes = new proxyHere[numNodes];
  int nProxyNodes;
  int i;

  // Note all home patches.
  for ( i = 0; i < numNodes; ++i )
  {
    proxyNodes[i] = No;
  }
  nProxyNodes=0;

  // Check all two-away neighbors.
  // This is really just one-away neighbors, since 
  // two-away always returns zero: RKB
  PatchID neighbors[PatchMap::MaxOneAway + PatchMap::MaxTwoAway];

  int myNode = patchMap->node(id);
  int numNeighbors = patchMap->downstreamNeighbors(id,neighbors);
  for ( i = 0; i < numNeighbors; ++i )
  {
    const int proxyNode = patchMap->node(neighbors[i]);
    if (proxyNode != myNode)
      if (proxyNodes[proxyNode] == No)
      {
	proxyNodes[proxyNode] = Yes;
	neighborNodes[nProxyNodes] = proxyNode;
	nProxyNodes++;
      }
  }

  delete [] proxyNodes;
  return nProxyNodes;
}

void LdbCoordinator::printLocalLdbReport(void)
{
  char outputBuf[255];
  char *curLoc;

  CkPrintf("%d:Patch report:\n",CkMyPe());
  
  curLoc = outputBuf;
  int i,j=0;
  for(i=0; i<patchMap->numPatches(); i++)
  {
    if (patchNAtoms[i] != -1)
    {
      curLoc += sprintf(curLoc,"%5d: %5d ",i,patchNAtoms[i]);
      j++;
    } 
    if (((j % 4) == 0) && j)
    {
      curLoc = outputBuf;
      CkPrintf("[%d]%s\n",CkMyPe(),outputBuf);
      j=0;
    }
  }

  CkPrintf("%d:Compute report:\n",CkMyPe());
  
  curLoc = outputBuf;
  j=0;
}

void LdbCoordinator::printRequiredProxies(PatchID id, FILE *fp)
{
  // Check all two-away neighbors.
  // This is really just one-away neighbors, since 
  // two-away always returns zero: RKB
  int neighborNodes[PatchMap::MaxOneAway + PatchMap::MaxTwoAway];
  const int nProxyNodes = requiredProxies(id,neighborNodes);

  fprintf(fp,"%4d ",nProxyNodes);

  for(int i=0;i<nProxyNodes;i++)
    fprintf(fp,"%4d ",neighborNodes[i]);
}

#include "LdbCoordinator.def.h"
