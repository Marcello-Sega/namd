
#if !defined(WIN32) || defined(__CYGWIN__)
#include <unistd.h>
#endif
#include <fcntl.h>

#include <charm++.h>

#include "NamdNborLB.h"
#include "NamdNborLB.def.h"
#include "Node.h"
#include "PatchMap.h"
#include "ComputeMap.h"
#include "LdbCoordinator.h"

void CreateNamdNborLB()
{
  // CkPrintf("[%d] creating NamdNborLB %d\n",CkMyPe(),loadbalancer);
  nborBaselb = CProxy_NamdNborLB::ckNew();
  // CkPrintf("[%d] created NamdNborLB %d\n",CkMyPe(),loadbalancer);
}

NamdNborLB::NamdNborLB()
{
  //  if (CkMyPe()==0)
  //   CkPrintf("[%d] NamdNborLB created\n",CkMyPe());
  processorArray = 0;
  patchArray = 0;
  computeArray = 0;
}

/*
NamdNborLB::~NamdNborLB()
{
  delete [] processorArray;
  delete [] patchArray;
  delete [] computeArray;
}
*/

CmiBool NamdNborLB::QueryBalanceNow(int _step)
{
  //  CkPrintf("[%d] Balancing on step %d\n",CkMyPe(),_step);
  if ( LdbCoordinator::Object()->takingLdbData ) {
    return CmiTrue;
  } else {
    return CmiFalse;
  }
}

NLBMigrateMsg* NamdNborLB::Strategy(NborBaseLB::LDStats* stats, int count)
{
#if CHARM_VERSION > 050403
  //  CkPrintf("LDB:[%d] All statistics received at %f, %f\n",
  //  CmiMyPe(), CmiTimer(),CmiWallTimer());
  int i,j;

  const int numProcessors = CmiNumPes();
  const int numPatches = PatchMap::Object()->numPatches();
  const int numComputes = ComputeMap::Object()->numComputes();
  const SimParameters* simParams = Node::Object()->simParameters;

  int nMoveableComputes = 0;
  for (i=0; i < count+1; i++) {
    LDStats &thisLDStats = ((i==count)?myStats:stats[i]);
    for (j=0; j < thisLDStats.n_objs; j++) {
      const LDObjData this_obj = thisLDStats.objData[j];
      if (this_obj.omID.id != 1) continue;
      if (this_obj.id.id[1] == -2) continue;
      if (this_obj.migratable)  nMoveableComputes++;
    }
  }
  CmiPrintf("%d nMoveableComputes: %d\n", CmiMyPe(), nMoveableComputes);

  // these sizes should never change
  processorArray = new processorInfo[numProcessors];
  patchArray = new patchInfo[numPatches];
//  if ( ! computeArray ) computeArray = new computeInfo[nMoveableComputes];
  computeArray = new computeInfo[nMoveableComputes];

  nMoveableComputes = buildData(stats,count);

  //CmiPrintf("AlgNbor begin on %d\n", CmiMyPe());
  AlgNbor(CkMyPe(), computeArray,patchArray,processorArray,
			nMoveableComputes, numPatches, numProcessors, count);
  //CmiPrintf("AlgNbor end on %d\n", CmiMyPe());
/*
  if (simParams->ldbStrategy == LDBSTRAT_REFINEONLY) {
    RefineOnly(computeArray,patchArray,processorArray,
                                nMoveableComputes, numPatches, numProcessors);
  } else if (simParams->ldbStrategy == LDBSTRAT_ALG7) {
    Alg7(computeArray,patchArray,processorArray,
                          nMoveableComputes, numPatches, numProcessors);
  } else if (simParams->ldbStrategy == LDBSTRAT_ALGROB) {
    if (step() == 0) {
      iout << iINFO << "Load balance cycle " << step()
        << " using RecBisection\n" << endi;
      AlgRecBisection(computeArray,patchArray,processorArray,
                            nMoveableComputes, numPatches, numProcessors);
    } else {
      iout << iINFO << "Load balance cycle " << step()
        << " using RefineOnly\n" << endi;
      RefineOnly(computeArray,patchArray,processorArray,
                                  nMoveableComputes, numPatches,
                                  numProcessors);
    }
  } else if (simParams->ldbStrategy == LDBSTRAT_OTHER) {
    if (step() == 0) {
      iout << iINFO << "Load balance cycle " << step()
        << " using Alg7\n" << endi;
      Alg7(computeArray,patchArray,processorArray,
                            nMoveableComputes, numPatches, numProcessors);
    } else {
      iout << iINFO << "Load balance cycle " << step()
        << " using RefineOnly\n" << endi;
      // To save the data to a file, uncomment the following lines -RKB
      //      if (step() == 1) {
      //	iout << iINFO << "Dumping data\n" << endi;
      //	dumpDataASCII("refinedata", numProcessors, numPatches,
      //		      nMoveableComputes);
      //      }
      RefineOnly(computeArray,patchArray,processorArray,
                                  nMoveableComputes, numPatches,
                                  numProcessors);
    }
  }
*/

  // For error checking:
  // Count up computes, to see if somebody doesn't have any computes
/*
  int* computeCount = new int[numProcessors];
  for(i=0; i<numProcessors; i++)
    computeCount[i]=0;
  for(i=0; i<nMoveableComputes; i++)
    computeCount[computeArray[i].processor]++;
  for(i=0; i<numProcessors; i++) {
    if (computeCount[i]==0)
      iout << iINFO <<"Warning: Processor " << i 
	   << " has NO moveable computes.\n" << endi;
  }
  delete [] computeCount;
*/
  
  CkVec<MigrateInfo *> migrateInfo;
  for(i=0;i<nMoveableComputes;i++) {
    if (computeArray[i].oldProcessor == CkMyPe())
    if (computeArray[i].processor != computeArray[i].oldProcessor) {
      // CkPrintf("[%d] Obj %d migrating from %d to %d\n",
      //          CkMyPe(),computeArray[i].handle.id.id[0],
      //       computeArray[i].processor,computeArray[i].oldProcessor);
      MigrateInfo *migrateMe = new MigrateInfo;
      migrateMe->obj = computeArray[i].handle;
      migrateMe->from_pe = computeArray[i].oldProcessor;
      migrateMe->to_pe = computeArray[i].processor;
      migrateInfo.insertAtEnd(migrateMe);
    }
  }
  
  int migrate_count=migrateInfo.length();
  CkPrintf("NamdNborLB [%d] migrating %d elements\n", CkMyPe(), migrate_count);
  NLBMigrateMsg* msg = new(&migrate_count,1) NLBMigrateMsg;
  msg->n_moves = migrate_count;
  for(i=0; i < migrate_count; i++) {
    MigrateInfo* item = migrateInfo[i];
    msg->moves[i] = *item;
    delete item;
    migrateInfo[i] = 0;
  }

  delete [] patchArray;
  delete [] computeArray;
  /*
  for(i=0; i<numProcessors; i++)
      delete [] processorArray[i].proxyUsage;
  */
  delete [] processorArray;

  return msg;
#else
  return NULL;
#endif
};


int NamdNborLB::buildData(NborBaseLB::LDStats* stats, int count)
{
#if CHARM_VERSION > 050403
  PatchMap* patchMap = PatchMap::Object();
  ComputeMap* computeMap = ComputeMap::Object();
  double bg_weight = 0.7;

  int i;
  for (i=0; i<CmiNumPes(); i++) {
    processorArray[i].load = 0.0;
    processorArray[i].backgroundLoad = 0.0;
    if (i == CmiMyPe()) {
      processorArray[i].Id = i;
    if (patchMap->numPatches() > 0)
      processorArray[i].backgroundLoad = myStats.bg_walltime*bg_weight;
    else 
      processorArray[i].backgroundLoad = myStats.bg_walltime;
      continue;
    }
    int peslot = NeighborIndex(i);
    if (peslot != -1) {
    processorArray[i].Id = i;
    if (patchMap->numPatches() > 0)
      processorArray[i].backgroundLoad = bg_weight * stats[peslot].bg_walltime;
    else 
      processorArray[i].backgroundLoad = stats[peslot].bg_walltime;
    }
    else 
      processorArray[i].Id = -2;
  }

  int nMoveableComputes=0;
  int nProxies = 0;		// total number of estimated proxies
  for (i=0; i<patchMap->numPatches(); i++) {
	patchArray[i].Id = i;
	patchArray[i].numAtoms = 0;
	patchArray[i].processor = patchMap->node(i);
	int neighborNodes[PatchMap::MaxOneAway + PatchMap::MaxTwoAway];
	const int numProxies = requiredProxies(i,neighborNodes);
        nProxies += numProxies;

	for (int k=0; k<numProxies; k++) {
	  if (NeighborIndex(neighborNodes[k]) != -1) {
	    processorArray[neighborNodes[k]].proxies.insert(&patchArray[i]);
	    patchArray[i].proxiesOn.insert(&processorArray[neighborNodes[k]]);
	  }
	}
  }
  for (i=0; i < count+1; i++) {
    int j;
    LDStats &thisLDStats = ((i==count)?myStats:stats[i]);
    for (j=0; j < thisLDStats.n_objs; j++) {
      const LDObjData this_obj = thisLDStats.objData[j];
      // filter out non-NAMD managed objects (like PME array)
      if (this_obj.omID.id != 1) continue;
      if (this_obj.id.id[1] == -2) { // Its a patch
/*
	const int pid = this_obj.id.id[0];
	int neighborNodes[PatchMap::MaxOneAway + PatchMap::MaxTwoAway];

	patchArray[pid].Id = pid;
	patchArray[pid].numAtoms = 0;
	patchArray[pid].processor = i;
	const int numProxies = requiredProxies(pid,neighborNodes);
        nProxies += numProxies;

	for (int k=0; k<numProxies; k++) {
	  processorArray[neighborNodes[k]].proxies.insert(&patchArray[pid]);
	  patchArray[pid].proxiesOn.insert(&processorArray[neighborNodes[k]]);
	}
*/
      } else if (this_obj.migratable) { // Its a compute
	const int cid = this_obj.id.id[0];
	const int p0 = computeMap->pid(cid,0);

	// For self-interactions, just return the same pid twice
	int p1;
	if (computeMap->numPids(cid) > 1)
	  p1 = computeMap->pid(cid,1);
	else p1 = p0;
	computeArray[nMoveableComputes].Id = cid;
	computeArray[nMoveableComputes].oldProcessor = thisLDStats.from_pe;
	computeArray[nMoveableComputes].processor = -1;
	computeArray[nMoveableComputes].patch1 = p0;
	computeArray[nMoveableComputes].patch2 = p1;
	computeArray[nMoveableComputes].handle = this_obj.handle;
	computeArray[nMoveableComputes].load = this_obj.wallTime;
	nMoveableComputes++;
      }
    }
  }
  return nMoveableComputes;
#else
  return 0;
#endif
}

// Figure out which proxies we will definitely create on other
// nodes, without regard for non-bonded computes.  This code is swiped
// from ProxyMgr, and changes there probable need to be propagated here.

int NamdNborLB::requiredProxies(PatchID id, int neighborNodes[])
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

  PatchMap* patchMap = PatchMap::Object();

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

