#include <iostream.h>
#include <stdlib.h>
#include "charm++.h"

#include "ProcessorPrivate.h"

#include "LdbCoordinator.decl.h"
#include "LdbCoordinator.h"
#include "Node.h"
#include "Namd.h"
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
#include "Alg7.h"
//#include "Alg0.h"
//#include "Alg1.h"
//#include "Alg4.h"
#include "packmsg.h"

#include "elements.h"
#include "ComputeMgr.decl.h"


class manheap;
class maxheap;

#define DEBUG_LEVEL 4

#define TIMER_FNC()   CmiWallTimer()

CmiHandler notifyIdleStart(void)
{
  LdbCoordinator::Object()->idleStart = TIMER_FNC();
  return 0;
}

CmiHandler notifyIdleEnd(void)
{
  const double idleEndTime = TIMER_FNC();
  LdbCoordinator::Object()->idleTime +=
    idleEndTime - LdbCoordinator::Object()->idleStart;
  LdbCoordinator::Object()->idleStart = -1;
  return 0;
}

LdbStatsMsg::LdbStatsMsg(void)
{
  nPatches = -1;
  nComputes = -1;
}

LdbStatsMsg::LdbStatsMsg(int np, int nc)
{
  nPatches = np;
  pid = new int[np];
  nAtoms = new int [np];
  nComputes = nc;
  cid = new int[nc];
  computeTime = new float[nc];
}

LdbStatsMsg::~LdbStatsMsg()
{
  delete [] pid;
  delete [] nAtoms;
  delete [] cid;
  delete [] computeTime;
}

PACK_MSG(LdbStatsMsg,
  PACK(proc);
  PACK(procLoad);
  PACK_AND_NEW_ARRAY(pid,nPatches);
  PACK_AND_NEW_ARRAY(nAtoms,nPatches);
  PACK_AND_NEW_ARRAY(cid,nComputes);
  PACK_AND_NEW_ARRAY(computeTime,nComputes);
)


LdbCoordinator::LdbCoordinator()
{
  if (CpvAccess(LdbCoordinator_instance) == NULL)
  {
    CpvAccess(LdbCoordinator_instance) = this;
  } else {
    iout << iFILE << iERROR << iPE 
	 << "LdbCoordinator instanced twice on same node!" << endi;
    CkExit();
  }

  ldbCycleNum = 1;
  nLocalComputes = nLocalPatches = 0;
  computeStartTime = computeTotalTime = (double *)NULL;
  patchNAtoms = (int *) NULL;
  sequencerThreads = (Sequencer **) NULL;
  if (CkMyPe() == 0)
    statsMsgs = new (LdbStatsMsg *[CkNumPes()]);
  else statsMsgs = NULL;
  ldbStatsFP = NULL;
  computeArray = NULL;
  patchArray = NULL;
  processorArray = NULL;

  // Register Converse timer routines
#ifndef NO_IDLE_COMPUTATION
  CsdSetNotifyIdle((CmiHandler)notifyIdleStart,(CmiHandler)notifyIdleEnd);
#endif
  idleTime = 0;
  nodesDone = 0;
}

LdbCoordinator::~LdbCoordinator(void)
{
  delete [] patchNAtoms;
  delete [] sequencerThreads;
  delete [] computeStartTime;
  delete [] computeTotalTime;
  if (CkMyPe() == 0)
  {
    delete [] statsMsgs;
    delete [] computeArray;
    delete [] patchArray;
    delete [] processorArray;
  }
  if (ldbStatsFP)
    fclose(ldbStatsFP);

}

void LdbCoordinator::initialize(PatchMap *pMap, ComputeMap *cMap, int reinit)
{
  const SimParameters *simParams = Node::Object()->simParameters;

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

  if ( ! reinit ) {
    delete [] sequencerThreads;  // Depends on delete NULL to do nothing
    sequencerThreads = new (Sequencer *[pMap->numPatches()]);
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
 
  delete [] computeStartTime;  // Depends on delete NULL to do nothing
  computeStartTime = new double[cMap->numComputes()];

  delete [] computeTotalTime;  // Depends on delete NULL to do nothing
  computeTotalTime = new double[cMap->numComputes()];

  nLocalComputes = 0;
  for(i=0;i<cMap->numComputes();i++)
  {
    if ( (cMap->node(i) == Node::Object()->myid())
	 && ( (cMap->type(i) == computeNonbondedPairType)
	      || (cMap->type(i) == computeNonbondedSelfType) ) )
    {
      computeStartTime[i] =
	computeTotalTime[i] = 0.;

      nLocalComputes++;
    } else {
      computeStartTime[i] = 
	computeTotalTime[i] = -1.;
    }
  }

  // Fixup to take care of the extra timestep at startup
  // This is pretty ugly here, but it makes the count correct
  if ((ldbCycleNum==1) || (ldbCycleNum == 2))
  {
    nLdbSteps = firstLdbStep;
  }
  else 
  {
    nLdbSteps = stepsPerLdbCycle;
  }

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
      processorArray = new processorInfo[Node::Object()->numNodes()];
  }
    

  if (nLocalPatches == 0 || nLocalComputes==0 )
  	checkAndGoToBarrier();

  // Start idle-time recording
  idleTime = 0;
#ifndef NO_IDLE_COMPUTATION
  CsdStartNotifyIdle();
#endif
  totalStartTime = TIMER_FNC();
}

void LdbCoordinator::patchLoad(PatchID id, int nAtoms, int /* timestep */)
{
  if (patchNAtoms[id] != -1)
  {
    patchNAtoms[id] = nAtoms;
    nPatchesReported++;
  } else {
    DebugM(10, "::patchLoad() Unexpected patch reporting in\n");
  }
  checkAndGoToBarrier();
}

void LdbCoordinator::startWork(ComputeID id, int /* timestep */ )
{
  if (Node::Object()->simParameters->ldbStrategy == LDBSTRAT_NONE)
    return;

  if (computeStartTime[id] == -1.)
  {
    DebugM(4, "::startWork() Unexpected compute("<<id<<") reporting in\n");
  } else if (computeStartTime[id] > 0.)
  {
    DebugM(4, "::startWork() Attempting to start already-running timer\n");
  }
  else
    computeStartTime[id] = TIMER_FNC();
}

void LdbCoordinator::endWork(ComputeID id, int /* timestep */)
{
  double endTime = TIMER_FNC();

  if (Node::Object()->simParameters->ldbStrategy == LDBSTRAT_NONE)
    return;

  if (computeStartTime[id] == -1.)
  {
    DebugM(4, "::endWork() Unexpected compute("<<id<<") reporting in\n");
  } else if (computeStartTime[id] == 0)
  {
    DebugM(4, "::endwork() Attempting to save non-running timer\n");
  } else  {
    computeTotalTime[id] += endTime-computeStartTime[id];
    computeStartTime[id] = 0.;
    nComputesReported++;
  }
  checkAndGoToBarrier();
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
  checkAndGoToBarrier();
  CthSuspend();
}

int LdbCoordinator::checkAndGoToBarrier(void)
{
  if ( (nPatchesReported > nPatchesExpected) 
       || (nComputesReported > nComputesExpected)
       || (controllerReported > controllerExpected) )
  {
    DebugM(3, "load balance barrier countdown: "
         << nPatchesReported << "/" << nPatchesExpected << " patches, "
         << nComputesReported << "/" << nComputesExpected << " computes, "
         << controllerReported << "/" << controllerExpected << " controller\n");
  }

  if (nPatchesReported == nPatchesExpected) 
  {
    DebugM(3, "All patches, " << nComputesReported << "/" <<
	nComputesExpected << " computes reported for load balance barrier.\n");
  }

  if ( (nPatchesReported == nPatchesExpected) 
       && (nComputesReported == nComputesExpected)
       && (controllerReported == controllerExpected) )
  {
    DebugM(3, "Load balance barrier reached.\n");
    LdbResumeMsg *msg = new LdbResumeMsg;
    CProxy_LdbCoordinator cl(thisgroup);
    cl.nodeDone(msg,0);
    return 1;
  }
  else return 0;
}

void LdbCoordinator::nodeDone(LdbResumeMsg *msg)
{
  nodesDone++;
  if (nodesDone==Node::Object()->numNodes())
  {
    nodesDone=0;
    CProxy_LdbCoordinator(thisgroup).sendStats(msg);
  }
  else delete msg;
}
void LdbCoordinator::sendStats(LdbResumeMsg *inMsg)
{
  delete inMsg;
/* REMOVE VERBOSE OUTPUT
  if(CkMyPe()==0)
  {
    CkPrintf("WallClock : %f  CPUTime : %f \n", CmiWallTimer()-Namd::cmiWallStart, 
	    CmiCpuTimer()-Namd::cmiCpuStart);
  }
*/
  // Turn off idle-time calculation
#ifndef NO_IDLE_COMPUTATION
  CsdStopNotifyIdle();
#endif
  totalTime = TIMER_FNC() - totalStartTime;

#ifndef NO_IDLE_COMPUTATION
  if (idleStart!= -1)
    iout << iPE << "WARNING: idle time still accumulating?\n" << endi;

#endif
  if (nLocalPatches > LDB_PATCHES)
  {
    char die_msg[255];
    sprintf(die_msg,
	    "%s(%d): Insufficient memory.  Increase LDB_PATCHES to %d",
	    __FILE__,__LINE__,nPatchesReported);
    NAMD_die(die_msg);
  }

  if (nLocalComputes > LDB_COMPUTES)
  {
    char die_msg[255];
    sprintf(die_msg,
	    "%s(%d): Insufficient memory.  Increase LDB_COMPUTES to %d",
	    __FILE__,__LINE__,nComputesReported);
    NAMD_die(die_msg);
  }

  LdbStatsMsg *msg = new LdbStatsMsg(nLocalPatches,nLocalComputes);
  
  if (msg == NULL)
    NAMD_die("LdbCoordinator::checkAndSendStats: Insufficient memory");

  msg->proc = Node::Object()->myid();
  msg->procLoad = totalTime - idleTime;

/*  REMOVE VERBOSE OUTPUT
  iout << iINFO << iPE << " Last " << nLdbSteps 
       << " steps: processor time = " << totalTime 
       << "  time per step = " << totalTime/nLdbSteps 
       << "\n" << endi;
#ifndef NO_IDLE_COMPUTATION
  CkPrintf("[%d] Processor idle time (this ldb cycle)=%5.1f%%\n",
	  msg->proc,100.*idleTime/totalTime);
#endif
*/

  int i;
  msg->nPatches = 0;

  for(i=0;i<patchMap->numPatches();i++)
  {
    if (patchNAtoms[i] != -1)
    {
      msg->pid[msg->nPatches]=i;
      msg->nAtoms[msg->nPatches]=patchNAtoms[i];
      msg->nPatches++;
    }
  }

  msg->nComputes = 0;
  for(i=0;i<computeMap->numComputes();i++)
  {
    if (computeStartTime[i] != -1.)
    {
      msg->cid[msg->nComputes]=i;
      msg->computeTime[msg->nComputes]=computeTotalTime[i];
      msg->nComputes++;
    }
  }
  for(i=0;i<msg->nComputes;i++)
    msg->procLoad -= msg->computeTime[i];
  
  CProxy_LdbCoordinator cl(thisgroup);
  cl.analyze(msg,0);
  return;
}

void LdbCoordinator::analyze(LdbStatsMsg *msg)
{
  if (Node::Object()->myid() != 0)
  {
    CkPrintf("Unexpected call to LdbCoordinator::analyze\n");
    return;
  }

  statsMsgs[msg->proc] = msg;
  nStatsMessagesReceived++;

  if (nStatsMessagesReceived==nStatsMessagesExpected)
  {
    processStatistics();
  }
}

void LdbCoordinator::processStatistics(void)
{
  //  CkPrintf("LDB: All statistics received at %f, %f\n",
  //  CmiTimer(),CmiWallTimer());

  const int numProcessors = Node::Object()->numNodes();
  const int numPatches = patchMap->numPatches();
  const SimParameters *simParams = Node::Object()->simParameters;

  const int nMoveableComputes = buildData();
  
  Rebalancer *rebalancer = 0;

  if (simParams->ldbStrategy == LDBSTRAT_REFINEONLY)
  {
    rebalancer = new RefineOnly(computeArray,patchArray,processorArray,
				nMoveableComputes, numPatches, numProcessors);
  } 
  else if (simParams->ldbStrategy == LDBSTRAT_ALG7)
  {
    rebalancer = new Alg7(computeArray,patchArray,processorArray,
			  nMoveableComputes, numPatches, numProcessors);
  }
  else if (simParams->ldbStrategy == LDBSTRAT_OTHER)
  {
    if (ldbCycleNum == 1)
    {
      iout << iINFO << "Load balance cycle " << ldbCycleNum
	<< " using Alg7\n" << endi;
      rebalancer = new Alg7(computeArray,patchArray,processorArray,
			    nMoveableComputes, numPatches, numProcessors);
    } else {
      iout << iINFO << "Load balance cycle " << ldbCycleNum
	<< " using RefineOnly\n" << endi;
      rebalancer = new RefineOnly(computeArray,patchArray,processorArray,
				  nMoveableComputes, numPatches,
				  numProcessors);
    }
  }


      
    
      
  delete rebalancer;

  // 0) Rebuild ComputeMap using computeMap->setNewNode()
  int i;
  for(i=0; i < nMoveableComputes; i++)
  {
    if ( (computeArray[i].processor != computeArray[i].oldProcessor)
	 && (computeArray[i].processor != -1) )
    {
      // CkPrintf("Assigning compute %d from %d to %d\n",
      //    i,computeArray[i].oldProcessor,computeArray[i].processor);
      computeMap->setNewNode(computeArray[i].Id,computeArray[i].processor);
      DebugM(2, "setting("<<computeArray[i].Id<<") newNode - curnode="
	<<computeMap->node(computeArray[i].Id) 
	<<" newnode="<<computeMap->newNode(computeArray[i].Id) << "\n" );
    } else {
      DebugM(2, "current setup - curnode="<<computeMap->node(computeArray[i].Id)
	<<" newnode="<<computeMap->newNode(computeArray[i].Id) << "\n" );
    }
  }
  /*
  if (computeMap->node(71) == 0) {
      iout << iPE << iERRORF << " moving computeID(71) by hand to node 2\n" << endi;
      computeMap->setNewNode(71, 2);
  }
  computeMap->setNewNode(72, (computeMap->node(72)+1)%3);
  computeMap->setNewNode(73, (computeMap->node(72)+1)%3);
  computeMap->setNewNode(74, (computeMap->node(72)+1)%3);
  computeMap->setNewNode(75, (computeMap->node(72)+2)%3);
  computeMap->setNewNode(76, (computeMap->node(72)+2)%3);
  computeMap->setNewNode(77, (computeMap->node(72)+2)%3);
  */

  // 1) Print out statistics in test format
  // printLdbReport(nMoveableComputes);


  // 2) delete messages
  cleanUpData();
  for (i=0; i < nStatsMessagesReceived; i++)
    delete statsMsgs[i];

  // computeMgr->updateComputes() call only on Node(0) i.e. right here
  // This will barrier for all Nodes - (i.e. Computes must be
  // here and with proxies before anyone can start up
  CProxy_ComputeMgr cm(CpvAccess(BOCclass_group).computeMgr);
  ComputeMgr *computeMgr = cm.ckLocalBranch();
  computeMgr->updateComputes(CProxy_LdbCoordinator::ckIdx_updateComputesReady(),thisgroup);
  //  CkPrintf("LDB: Done processing statistics at %f, %f\n",
  //	  CmiTimer(),CmiWallTimer());
}

void LdbCoordinator::updateComputesReady() {
  DebugM(3,"updateComputesReady()\n");

  LdbResumeMsg *sendmsg = new LdbResumeMsg;
  CProxy_LdbCoordinator(thisgroup).resume(sendmsg);
  CkStartQD(CProxy_LdbCoordinator::ckIdx_resumeReady((CkQdMsg*)0),&thishandle);
}

void LdbCoordinator::resume(LdbResumeMsg *msg)
{
  DebugM(3,"resume()\n");
  //  printLocalLdbReport();

  ldbCycleNum++;
  initialize(PatchMap::Object(),ComputeMap::Object(),1);
  delete msg;
}

void LdbCoordinator::resumeReady(CkQdMsg *msg) {
  DebugM(3,"resumeReady()\n");
  delete msg;

  Namd::startTimer();
  LdbResumeMsg *sendmsg = new LdbResumeMsg;
  CProxy_LdbCoordinator(thisgroup).resume2(sendmsg);
}

void LdbCoordinator::resume2(LdbResumeMsg *msg)
{
  DebugM(3,"resume2()\n");
  awakenSequencers();
  delete msg;
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
      //      CkPrintf("Awakening thread %d\n",i);
      sequencerThreads[i]->awaken();
    }
    sequencerThreads[i]= NULL;
  }
}

int LdbCoordinator::buildData(void)
{
  int i;
  for (i=0; i<nStatsMessagesReceived; i++)
  {
    const LdbStatsMsg *msg = statsMsgs[i];
    processorArray[i].Id = msg->proc;
    processorArray[i].backgroundLoad = statsMsgs[i]->procLoad;
    processorArray[i].proxies = new Set();
  }

  int nMoveableComputes=0;
  for (i=0; i < nStatsMessagesReceived; i++)
  {
    const LdbStatsMsg *msg = statsMsgs[i];
    int j;
    for (j=0; j<msg->nPatches; j++)
    {
      const int pid = msg->pid[j];
      int neighborNodes[PatchMap::MaxOneAway + PatchMap::MaxTwoAway];

      patchArray[pid].Id = pid;
      patchArray[pid].numAtoms = msg->nAtoms[j];
      patchArray[pid].processor = patchMap->node(pid);
      const int numProxies = requiredProxies(pid,neighborNodes);
      patchArray[pid].proxiesOn = new Set();

      for (int k=0; k<numProxies; k++)
      {
	processorArray[neighborNodes[k]].proxies->insert(&patchArray[pid]);
	patchArray[pid].proxiesOn->insert(&processorArray[neighborNodes[k]]);
      }
    }

    for (j=0; j<statsMsgs[i]->nComputes; j++)
    {
      const int p0 = computeMap->pid(msg->cid[j],0);

      // For self-interactions, just return the same pid twice
      int p1;
      if (computeMap->numPids(msg->cid[j]) > 1)
	p1 = computeMap->pid(msg->cid[j],1);
      else 
	p1 = p0;

      const int cid = msg->cid[j];
      computeArray[nMoveableComputes].Id = cid;
      computeArray[nMoveableComputes].oldProcessor = msg->proc;
      computeArray[nMoveableComputes].patch1 = p0;
      computeArray[nMoveableComputes].patch2 = p1;
      computeArray[nMoveableComputes].load = msg->computeTime[j];
      nMoveableComputes++;
    }
/* REMOVE VERBOSE OUTPUT
    CkPrintf("PE %d nComputes = %d\n",msg->proc,j);
*/
    
  }
  return nMoveableComputes;
}

void LdbCoordinator::cleanUpData(void)
{
  const int numPatches = patchMap->numPatches();
  
  int i;
  for (i=0; i<nStatsMessagesReceived; i++)
    delete processorArray[i].proxies;

  for (i=0; i < numPatches; i++)
    delete patchArray[i].proxiesOn;

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
  for(i=0; i<computeMap->numComputes(); i++)
  {
    if (computeTotalTime[i] != -1)
    {
      curLoc += sprintf(curLoc,"%5d: %4f ",i,computeTotalTime[i]);
      j++;
    } 
    if (((j % 4) == 0) && j)
    {
      curLoc = outputBuf;
      CkPrintf("[%d]%s\n",CkMyPe(),outputBuf);
      j=0;
    }
  }
    
}

void LdbCoordinator::printLdbReport(const int nMoveableComputes)
{
  if (ldbStatsFP == NULL)
  {
    ldbStatsFP = fopen("ldbreport.dat","w");
  }

  const int nProcs = Node::Object()->numNodes();
  const int nPatches = patchMap->numPatches();

  int i,j;

  fprintf(ldbStatsFP,"*** Load balancer report ***\n");

  fprintf(ldbStatsFP,"%4d %4d %4d\n",nProcs,nPatches,nMoveableComputes);

  // Print out processor background load
  for (i=0;i<nProcs;i++)
    fprintf(ldbStatsFP,"%8.3f ",statsMsgs[i]->procLoad);
  fprintf(ldbStatsFP,"\n\n");

  // Print out info for each patch
  for (i=0;i<nProcs;i++)
  {
    const LdbStatsMsg *msg = statsMsgs[i];
    for ( j=0; j < msg->nPatches; j++)
    {
      fprintf(ldbStatsFP,"%4d %4d %4d ",
	      msg->pid[j],msg->nAtoms[j],msg->proc);
      printRequiredProxies(msg->pid[j],ldbStatsFP);
      fprintf(ldbStatsFP,"\n");
    }
  }
  fprintf(ldbStatsFP,"\n");
  
  // Print out info for each compute
  int numComputesPrinted = 0;
  for (i=0;i<nProcs;i++)
  {
    const LdbStatsMsg *msg = statsMsgs[i];
    for ( j=0; j < msg->nComputes; j++)
    {
      const int p0 = computeMap->pid(msg->cid[j],0);

      // For self-interactions, just return the same pid twice
      int p1;
      if (computeMap->numPids(msg->cid[j]) > 1)
	p1 = computeMap->pid(msg->cid[j],1);
      else 
	p1 = p0;

      fprintf(ldbStatsFP,"%4d %4d %4d %4d %8.3f\n",
	      msg->cid[j],msg->proc,p0,p1,msg->computeTime[j]);
      numComputesPrinted++;
    }
  }
  if (numComputesPrinted != nMoveableComputes)
    CkPrintf("LDB stats missing %d %d \n",
	    numComputesPrinted,nMoveableComputes);
  fprintf(ldbStatsFP,"\n");
  fprintf(ldbStatsFP,"*** Load balancer report complete ***\n");
  fflush(ldbStatsFP);
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
