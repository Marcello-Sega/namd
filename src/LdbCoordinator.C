#include <iostream.h>
#include <stdlib.h>
#include "ckdefs.h"
#include "chare.h"
#include "c++interface.h"

#include "LdbCoordinator.top.h"
#include "LdbCoordinator.h"
#include "Node.h"
#include "SimParameters.h"
#include "PatchMap.h"
#include "ComputeMap.h"
#define DEBUGM
#include "Debug.h"
#include "Sequencer.h"
#include "RefineOnly.h"
#include "ComputeMgr.h"
#include "Alg7.h"
//#include "Alg0.h"
//#include "Alg1.h"
//#include "Alg4.h"

#include "elements.h"
#include "RefineOnly.h"


class manheap;
class maxheap;

#define DEBUG_LEVEL 3
#define TIMER_FNC()   CmiTimer()

// static initialization
LdbCoordinator *LdbCoordinator::_instance = 0;

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

LdbCoordinator::LdbCoordinator(InitMsg *msg)
{
  if (_instance == NULL)
  {
    _instance = this;
  } else {
    iout << iFILE << iERROR << iPE 
	 << "LdbCoordinator instanced twice on same node!" << endi;
    CharmExit();
  }

  first_ldbcycle = TRUE;
  nLocalComputes = nLocalPatches = 0;
  computeStartTime = computeTotalTime = (double *)NULL;
  patchNAtoms = (int *) NULL;
  sequencerThreads = (Sequencer **) NULL;
  if (CMyPe() == 0)
    statsMsgs = new (LdbStatsMsg *[CNumPes()]);
  else statsMsgs = NULL;
  ldbStatsFP = NULL;
  computeArray = NULL;
  patchArray = NULL;
  processorArray = NULL;

  delete msg;

  // Register Converse timer routines
  //CsdSetNotifyIdle((CmiHandler)notifyIdleStart,(CmiHandler)notifyIdleEnd);
  idleTime = 0;
}

LdbCoordinator::~LdbCoordinator(void)
{
  delete [] patchNAtoms;
  delete [] sequencerThreads;
  delete [] computeStartTime;
  delete [] computeTotalTime;
  if (CMyPe() == 0)
  {
    delete [] statsMsgs;
    delete [] computeArray;
    delete [] patchArray;
    delete [] processorArray;
  }
  if (ldbStatsFP)
    fclose(ldbStatsFP);

}

void LdbCoordinator::initialize(PatchMap *pMap, ComputeMap *cMap)
{
  //  DebugM(10,"stepsPerLdbCycle initialized\n");
  stepsPerLdbCycle = (Node::Object()->simParameters)->ldbStepsPerCycle;

  computeMap = cMap;
  patchMap = pMap;

  // Set the number of received messages correctly for node 0

  nStatsMessagesExpected = Node::Object()->numNodes();
  nStatsMessagesReceived = 0;

  delete [] patchNAtoms;  // Depends on delete NULL to do nothing
  patchNAtoms = new int[pMap->numPatches()];

  delete [] sequencerThreads;  // Depends on delete NULL to do nothing
  sequencerThreads = new (Sequencer *[pMap->numPatches()]);

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
    sequencerThreads[i]=NULL;
  }
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

  nPatchesReported = 0;
  nPatchesExpected = nLocalPatches;

  nComputesReported = 0;
  nComputesExpected = nLocalComputes * stepsPerLdbCycle;
  // Fixup to take care of the extra timestep at startup
  // This is pretty ugly here, but it makes the count correct
  if (first_ldbcycle)
  {
    nComputesExpected += nLocalComputes;
    first_ldbcycle = FALSE;
  }

  if (CMyPe() == 0)
  {
    if (computeArray == NULL)
      computeArray = new computeInfo[computeMap->numComputes()];
    if (patchArray == NULL)
      patchArray = new patchInfo[patchMap->numPatches()];
    if (processorArray == NULL)
      processorArray = new processorInfo[Node::Object()->numNodes()];
  }
    

  // Start idle-time recording
  idleTime = 0;
  //  CsdStartNotifyIdle();
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
  checkAndSendStats();
}

void LdbCoordinator::startWork(ComputeID id, int /* timestep */ )
{
  if (Node::Object()->simParameters->ldbStrategy == LDBSTRAT_NONE)
    return;

  if (computeStartTime[id] == -1.)
  {
    DebugM(4, "::startWork() Unexpected compute reporting in\n");
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
    DebugM(4, "::endWork() Unexpected compute reporting in\n");
  } else if (computeStartTime[id] == 0)
  {
    DebugM(4, "::endwork() Attempting to save non-running timer\n");
  } else  {
    computeTotalTime[id] += endTime-computeStartTime[id];
    computeStartTime[id] = 0.;
    nComputesReported++;
  }
  checkAndSendStats();
}

void LdbCoordinator::rebalance(Sequencer *seq, PatchID pid)
{
  if (Node::Object()->simParameters->ldbStrategy == LDBSTRAT_NONE)
    return;

  sequencerThreads[pid] = seq;
  seq->suspend();
}

int LdbCoordinator::checkAndSendStats(void)
{
  if ( (nPatchesReported == nPatchesExpected) 
       && (nComputesReported == nComputesExpected) )
  {
    // Turn off idle-time calculation
    //CsdStopNotifyIdle();
    totalTime = TIMER_FNC() - totalStartTime;
    if (idleStart!= -1)
      iout << iPE << "WARNING: idle time still accumulating?\n" << endi;

    // Here, all the data gets sent to Node 0
    // For now, just send a dummy message to Node 0.

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

    LdbStatsMsg *msg = new (MsgIndex(LdbStatsMsg)) LdbStatsMsg;

    if (msg == NULL)
      NAMD_die("LdbCoordinator::checkAndSendStats: Insufficient memory");

    msg->proc = Node::Object()->myid();
    msg->procLoad = totalTime - idleTime;
    CPrintf("[%d] Processor idle time=%5.1f%%\n",
	    msg->proc,100.*idleTime/totalTime);
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
    for(i=0;i<nLocalComputes;i++)
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

    CSendMsgBranch(LdbCoordinator, analyze, msg, thisgroup,0);
    return 1;
  }
  else return 0;
}

void LdbCoordinator::analyze(LdbStatsMsg *msg)
{
  CPrintf("Node %d receiving LDB results\n",CMyPe());

  if (Node::Object()->myid() != 0)
  {
    CPrintf("Unexpected call to LdbCoordinator::analyze\n");
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
  CPrintf("All statistics received\n");

  const int numProcessors = Node::Object()->numNodes();
  const int numPatches = patchMap->numPatches();

  const int nMoveableComputes = buildData();
  
  char *algChoice = "Alg7"; 

  if (strcmp(algChoice,"RefineOnly") == 0)
  {
    new RefineOnly(computeArray,patchArray,processorArray,
		   nMoveableComputes, numPatches, numProcessors);
  } 
  else if (strcmp(algChoice,"Alg7") == 0)
  {
    new Alg7(computeArray,patchArray,processorArray,
	     nMoveableComputes, numPatches, numProcessors);
  }
//   else if (algChoice == "Alg1")
//     Alg1::Alg1(computeArray,patchArray,processorArray,
//       numComputes, numPatches, numProcessors);
//   else if (algChoice == "Alg4")
//     Alg4::Alg4(computeArray,patchArray,processorArray,
//       numComputes, numPatches, numProcessors);
//   else if (algChoice == "Alg1")
//     Alg7::Alg7(computeArray,patchArray,processorArray,
//       numComputes, numPatches, numProcessors);

  // 0) Rebuild ComputeMap using computeMap->setNewNode()
  int i;
  for(i=0; i < nMoveableComputes; i++)
  {
    if ( (computeArray[i].processor != computeArray[i].oldProcessor)
	 && (computeArray[i].processor != -1) )
    {
      // CPrintf("Assigning compute %d from %d to %d\n",
      //    i,computeArray[i].oldProcessor,computeArray[i].processor);
      computeMap->setNewNode(computeArray[i].Id,computeArray[i].processor);
    }
  }

  // 1) Print out statistics in test format
  printLdbReport(nMoveableComputes);

  // 2) delete messages
  cleanUpData();
  for (i=0; i < nStatsMessagesReceived; i++)
    delete statsMsgs[i];

  // computeMgr->updateComputes() call only on Node(0) i.e. right here
  // This will barrier for all Nodes - (i.e. Computes must be
  // here and with proxies before anyone can start up

  ComputeMgr *computeMgr = CLocalBranch(ComputeMgr, group.computeMgr);
  computeMgr->updateComputes(GetEntryPtr(LdbCoordinator,updateComputesReady),thisgroup);
}

void LdbCoordinator::updateComputesReady(DoneMsg *msg) {
  delete msg;

  LdbResumeMsg *sendmsg = new (MsgIndex(LdbResumeMsg)) LdbResumeMsg;
  CBroadcastMsgBranch(LdbCoordinator, resume, sendmsg, thisgroup);
}

void LdbCoordinator::resume(LdbResumeMsg *msg)
{
  //  printLocalLdbReport();

  awakenSequencers();
  initialize(PatchMap::Object(),ComputeMap::Object());
  delete msg;
}

void LdbCoordinator::awakenSequencers()
{
  for(int i=0; i < patchMap->numPatches(); i++)
  {
    if (sequencerThreads[i])
    {
      //      CPrintf("Awakening thread %d\n",i);
      sequencerThreads[i]->awaken();
    }
    sequencerThreads[i]= NULL;
  }
}

int LdbCoordinator::buildData(void)
{
  const int numProcessors = Node::Object()->numNodes();
  const int numPatches = patchMap->numPatches();
  
  int i;
  for (i=0; i<nStatsMessagesReceived; i++)
  {
    const LdbStatsMsg *msg = statsMsgs[i];
    processorArray[i].Id = msg->proc;
    processorArray[i].backgroundLoad = 0;
                             // should be = statsMsgs[i]->procLoad;
    processorArray[i].proxies = new Set();
  }

  int nMoveableComputes=0;
  for (i=0; i < nStatsMessagesReceived; i++)
  {
    const LdbStatsMsg *msg = statsMsgs[i];
    for (int j=0; j<msg->nPatches; j++)
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
      computeArray[cid].Id = cid;
      computeArray[cid].oldProcessor = msg->proc;
      computeArray[cid].patch1 = p0;
      computeArray[cid].patch2 = p1;
      computeArray[cid].load = msg->computeTime[j];
      nMoveableComputes++;
    }
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
  int numPatches = patchMap->numPatches();
  enum proxyHere { No, Yes };
  int numNodes = CNumPes();
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
  int numNeighbors = patchMap->oneAwayNeighbors(id,neighbors);
  numNeighbors += patchMap->twoAwayNeighbors(id,neighbors+numNeighbors);
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

  CPrintf("%d:Patch report:\n",CMyPe());
  
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
      CPrintf("[%d]%s\n",CMyPe(),outputBuf);
      j=0;
    }
  }

  CPrintf("%d:Compute report:\n",CMyPe());
  
  curLoc = outputBuf;
  j=0;
  for(i=0; i<computeMap->numComputes(); i++)
  {
    if (computeTotalTime[i] != -1)
    {
      curLoc += sprintf(curLoc,"%5d: %4lf ",i,computeTotalTime[i]);
      j++;
    } 
    if (((j % 4) == 0) && j)
    {
      curLoc = outputBuf;
      CPrintf("[%d]%s\n",CMyPe(),outputBuf);
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
    CPrintf("LDB stats missing %d %d \n",
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

#include "LdbCoordinator.bot.h"
