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

  nLocalComputes = nLocalPatches = 0;
  computeStartTime = computeTotalTime = (double *)NULL;
  patchNAtoms = (int *) NULL;
  sequencerThreads = (Sequencer **) NULL;
  delete msg;
}

LdbCoordinator::~LdbCoordinator(void)
{
  delete [] patchNAtoms;
  delete [] sequencerThreads;
  delete [] computeStartTime;
  delete [] computeTotalTime;
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
  sequencerThreads[pid] = seq;
  seq->suspend();
}

int LdbCoordinator::checkAndSendStats(void)
{
  if ( (nPatchesReported >= nPatchesExpected) 
       && (nComputesReported >= nComputesExpected) )
  {
    CPrintf("nPatchesReported = %d, nComputesReported = %d\n",
	    nPatchesReported-nPatchesExpected, 
	    nComputesReported-nComputesExpected);
    // Here, all the data gets sent to Node 0

    // For now, no send to Node 0.  Instead just call resume() 
    // with a dummy message to avoid thread problems.

    LdbStatsMsg *msg = new (MsgIndex(LdbStatsMsg)) LdbStatsMsg;
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

  nStatsMessagesReceived++;
  if (nStatsMessagesReceived==nStatsMessagesExpected)
  {
    CPrintf("Node 0 LDB resuming other processors\n",CMyPe());
    LdbResumeMsg *sendmsg = new (MsgIndex(LdbResumeMsg)) LdbResumeMsg;
    CBroadcastMsgBranch(LdbCoordinator, resume, sendmsg, thisgroup);
  }
  delete msg;
}

void LdbCoordinator::resume(LdbResumeMsg *msg)
{
  printLocalLdbReport();

  awakenSequencers();
  initialize(PatchMap::Object(),ComputeMap::Object());
  delete msg;
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

#include "LdbCoordinator.bot.h"
