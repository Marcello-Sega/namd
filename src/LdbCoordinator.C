#include "ckdefs.h"
#include "chare.h"
#include "c++interface.h"

#include "LdbCoordinator.top.h"
#include "LdbCoordinator.h"
#include "Node.h"
#include "SimParameters.h"
#include "PatchMap.h"
#include "ComputeMap.h"
#include "Debug.h"
#include "Sequencer.h"

#define DEBUG_LEVEL 3
// static initialization
LdbCoordinator *LdbCoordinator::_instance = 0;

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

  delete msg;
}

LdbCoordinator::~LdbCoordinator(void)
{
  delete [] patchNAtoms;
  delete [] sequencerThreads;
  delete [] computeStartTime;
  delete [] computeTotalTime;
  if (CMyPe() == 0)
    delete [] statsMsgs;
}

void LdbCoordinator::initialize(PatchMap *pMap, ComputeMap *cMap)
{
  DebugM(10,"stepsPerLdbCycle initialized\n");
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
}

void LdbCoordinator::patchLoad(PatchID id, int nAtoms, int timestep)
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

void LdbCoordinator::startWork(ComputeID id, int timestep)
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
    computeStartTime[id] = CmiTimer();
}

void LdbCoordinator::endWork(ComputeID id, int timestep)
{
  double endTime = CmiTimer();

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
    msg->procLoad = 0;

    msg->nPatches = 0;

    int i;
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
    CPrintf("All statistics received\n");
    // 1) Print out statistics in test format
    // 2) Delete messages
    // 3) Resume operation with a message
    
    // 1) Print out statistics in test format
    printLdbReport();

    // 2) delete messages
    for (int i=0; i < nStatsMessagesReceived; i++)
      delete statsMsgs[i];

    // 3) Resume operation with a message
    CPrintf("Node 0 LDB resuming other processors\n",CMyPe());
    LdbResumeMsg *sendmsg = new (MsgIndex(LdbResumeMsg)) LdbResumeMsg;
    CBroadcastMsgBranch(LdbCoordinator, resume, sendmsg, thisgroup);
  }
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

void LdbCoordinator::printLdbReport(void)
{
  const int nProcs = Node::Object()->numNodes();
  const int nPatches = patchMap->numPatches();
  const int nComputes = computeMap->numComputes();

  int i,j;

  CPrintf("*** Load balancer report ***\n");

  CPrintf("%4d %4d %4d\n",nProcs,nPatches,nComputes);

  // Print out processor background load
  for (i=0;i<nProcs;i++)
    CPrintf("%4d ",statsMsgs[i]->procLoad);
  CPrintf("\n\n");

  // Print out info for each patch
  for (i=0;i<nProcs;i++)
  {
    const LdbStatsMsg *msg = statsMsgs[i];
    for ( j=0; j < msg->nPatches; j++)
      CPrintf("%4d %4d %4d %4d\n",msg->pid[j],msg->nAtoms[j],msg->proc,0);
  }
  CPrintf("\n");
  
  // Print out info for each compute
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

      CPrintf("%4d %4d %4d %4d %8.3f\n",
	      msg->cid[j],msg->proc,p0,p1,msg->computeTime[j]);
    }
  }
  CPrintf("\n");
  CPrintf("*** Load balancer report complete ***\n");
}
#include "LdbCoordinator.bot.h"
