
#ifndef WIN32
#include <unistd.h>
#endif
#include <fcntl.h>

#include <charm++.h>

#include <CkLists.h>
#include "NamdCentLB.h"
#include "NamdCentLB.def.h"
#include "Node.h"
#include "PatchMap.h"
#include "ComputeMap.h"

void CreateNamdCentLB()
{
  // CkPrintf("[%d] creating NamdCentLB %d\n",CkMyPe(),loadbalancer);
  loadbalancer = CProxy_NamdCentLB::ckNew();
  // CkPrintf("[%d] created NamdCentLB %d\n",CkMyPe(),loadbalancer);
}

NamdCentLB::NamdCentLB()
{
  //  if (CkMyPe()==0)
  //   CkPrintf("[%d] NamdCentLB created\n",CkMyPe());
}

CmiBool NamdCentLB::QueryBalanceNow(int _step)
{
  //  CkPrintf("[%d] Balancing on step %d\n",CkMyPe(),_step);
  return CmiTrue;
}

CLBMigrateMsg* NamdCentLB::Strategy(CentralLB::LDStats* stats, int count)
{
  //  CkPrintf("LDB: All statistics received at %f, %f\n",
  //  CmiTimer(),CmiWallTimer());

  const int numProcessors = count;
  const int numPatches = PatchMap::Object()->numPatches();
  const int numComputes = ComputeMap::Object()->numComputes();
  const SimParameters* simParams = Node::Object()->simParameters;

  processorArray = new processorInfo[numProcessors];
  patchArray = new patchInfo[numPatches];
  computeArray = new computeInfo[numComputes];

  const int nMoveableComputes = buildData(stats,count);
  // gzheng debug
  // dumpData("data", numProcessors, numPatches, nMoveableComputes);
  // loadData("data", numProcessors, numPatches, nMoveableComputes);
  // end of debug section

  Rebalancer* rebalancer = 0;

  if (simParams->ldbStrategy == LDBSTRAT_REFINEONLY) {
    rebalancer = new RefineOnly(computeArray,patchArray,processorArray,
                                nMoveableComputes, numPatches, numProcessors);
  } else if (simParams->ldbStrategy == LDBSTRAT_ALG7) {
    rebalancer = new Alg7(computeArray,patchArray,processorArray,
                          nMoveableComputes, numPatches, numProcessors);
  } else if (simParams->ldbStrategy == LDBSTRAT_OTHER) {
    if (step() == 0) {
      iout << iINFO << "Load balance cycle " << step()
        << " using Alg7\n" << endi;
      rebalancer = new Alg7(computeArray,patchArray,processorArray,
                            nMoveableComputes, numPatches, numProcessors);
    } else {
      iout << iINFO << "Load balance cycle " << step()
        << " using RefineOnly\n" << endi;
      rebalancer = new RefineOnly(computeArray,patchArray,processorArray,
                                  nMoveableComputes, numPatches,
                                  numProcessors);
    }
  }

  CkVector migrateInfo;
  int i;
  for(i=0;i<nMoveableComputes;i++) {
    if (computeArray[i].processor != computeArray[i].oldProcessor) {
      //      CkPrintf("[%d] Obj %d migrating from %d to %d\n",
      //               CkMyPe(),computeArray[i].handle.id.id[0],
      //	       computeArray[i].processor,computeArray[i].oldProcessor);
      MigrateInfo *migrateMe = new MigrateInfo;
      migrateMe->obj = computeArray[i].handle;
      migrateMe->from_pe = computeArray[i].oldProcessor;
      migrateMe->to_pe = computeArray[i].processor;
      migrateInfo.push_back((void*)migrateMe);
    }
  }
  
  delete [] processorArray;
  delete [] patchArray;
  delete [] computeArray;

  int migrate_count=migrateInfo.size();
  // CkPrintf("NamdCentLB migrating %d elements\n",migrate_count);
  CLBMigrateMsg* msg = new(&migrate_count,1) CLBMigrateMsg;
  msg->n_moves = migrate_count;
  for(i=0; i < migrate_count; i++) {
    MigrateInfo* item = (MigrateInfo*) migrateInfo[i];
    msg->moves[i] = *item;
    delete item;
    migrateInfo[i] = 0;
  }
  return msg;
};

#ifndef WIN32

void NamdCentLB::dumpData(char *file, int numProcessors, int numPatches, int numComputes)
{
  int fd = open(file, O_WRONLY|O_CREAT|O_TRUNC, 0644);
  if (fd == -1){
     perror("dumpLDStats");
     return;
  }
  CkPrintf("***** DUMP data to file: %s ***** \n", file);
  write(fd, &numProcessors, sizeof(int));
  write(fd, &numPatches, sizeof(int));
  write(fd, &numComputes, sizeof(int));
  write(fd, processorArray, sizeof(processorInfo)*numProcessors);
  write(fd, patchArray, sizeof(patchInfo)*numPatches);
  write(fd, computeArray, sizeof(computeInfo)*numComputes);
  int i;
  // dump patchSet
  for (i=0; i< numProcessors; i++)
  {
      int num = processorArray[i].proxies->numElements();
      write(fd, &num, sizeof(int));
//CkPrintf("**** Proc:%d num:%d \n", i, num);
      Iterator nextProxy;
//      nextProxy.id = 0;
      patchInfo *p = (patchInfo *)processorArray[i].proxies->iterator((Iterator *)&nextProxy);
      while (p) 
      {
          write(fd, &p->Id, sizeof(int));
          p = (patchInfo *)processorArray[i].proxies->next((Iterator*)&nextProxy);
      }
  }
  // dump proxiesOn
  for (i=0; i<numPatches; i++)
  {
      int num = patchArray[i].proxiesOn->numElements();
      write(fd, &num, sizeof(int));
// CkPrintf("**** Patch:%d num:%d \n", i, num);
      Iterator nextProc;
//      nextProc.id = 0;
      processorInfo *p = (processorInfo *)patchArray[i].proxiesOn->iterator((Iterator *)&nextProc);
      while (p) 
      {
          write(fd, &p->Id, sizeof(int));
          p = (processorInfo *)patchArray[i].proxiesOn->next((Iterator*)&nextProc);
      }
  }

  close(fd);
}

void NamdCentLB::loadData(char *file, int &numProcessors, int &numPatches, int &numComputes)
{
  int fd = open(file, O_RDONLY);
  if (fd == -1){
     perror("loadData");
     return;
  }

  read(fd, &numProcessors, sizeof(int));
  read(fd, &numPatches, sizeof(int));
  read(fd, &numComputes, sizeof(int));

  printf("numProcs: %d numPatches: %d numComputes: %d\n", numProcessors,numPatches, numComputes);
  processorInfo *processorArray = new processorInfo[numProcessors];
  computeInfo *computeArray = new computeInfo[numComputes];
  patchInfo *patchArray = new patchInfo [numPatches];

  read(fd, processorArray, sizeof(processorInfo)*numProcessors);
  read(fd, patchArray, sizeof(patchInfo)*numPatches);
  read(fd, computeArray, sizeof(computeInfo)*numComputes);

  int i;
  for (i=0; i<numProcessors; i++)
  {
      int num;
      read(fd, &num, sizeof(int));
      printf("proc: %d proxies:%d \n", i, num);
      processorArray[i].proxies = new Set();
      for (int j=0; j<num; j++) {
          int id;
          read(fd, &id, sizeof(int));
	  printf("%d ", id);
          processorArray[i].proxies->insert(&patchArray[id]);
      }
      printf("\n");
  }
  for (i=0; i<numPatches; i++)
  {
      int num;
      read(fd, &num, sizeof(int));
//      printf("patch: %d proxiesOn:%d \n", i, num);
      patchArray[i].proxiesOn = new Set();
      for (int j=0; j<num; j++) {
          int id;
          read(fd, &id, sizeof(int));
          patchArray[i].proxiesOn->insert(&processorArray[id]);
      }
  }

  close(fd);
}

#endif

int NamdCentLB::buildData(CentralLB::LDStats* stats, int count)
{
  PatchMap* patchMap = PatchMap::Object();
  ComputeMap* computeMap = ComputeMap::Object();
  double bg_weight = 1.0;

  int i;
  for (i=0; i<count; i++) {
    processorArray[i].Id = i;
    processorArray[i].backgroundLoad = bg_weight * stats[i].bg_walltime;
    processorArray[i].proxies = new Set();
  }

  int nMoveableComputes=0;
  for (i=0; i < count; i++) {
    int j;
    for (j=0; j < stats[i].n_objs; j++) {
      const LDObjData this_obj = stats[i].objData[j];
      if (this_obj.id.id[1] == -2) { // Its a patch
	const int pid = this_obj.id.id[0];
	int neighborNodes[PatchMap::MaxOneAway + PatchMap::MaxTwoAway];

	patchArray[pid].Id = pid;
	patchArray[pid].numAtoms = 0;
	patchArray[pid].processor = i;
	const int numProxies = requiredProxies(pid,neighborNodes);
	patchArray[pid].proxiesOn = new Set();

	for (int k=0; k<numProxies; k++) {
	  processorArray[neighborNodes[k]].proxies->insert(&patchArray[pid]);
	  patchArray[pid].proxiesOn->insert(&processorArray[neighborNodes[k]]);
	}
      } else if (this_obj.migratable) { // Its a compute
	const int cid = this_obj.id.id[0];
	const int p0 = computeMap->pid(cid,0);

	// For self-interactions, just return the same pid twice
	int p1;
	if (computeMap->numPids(cid) > 1)
	  p1 = computeMap->pid(cid,1);
	else p1 = p0;
	computeArray[nMoveableComputes].Id = cid;
	computeArray[nMoveableComputes].oldProcessor = i;
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
}

// Figure out which proxies we will definitely create on other
// nodes, without regard for non-bonded computes.  This code is swiped
// from ProxyMgr, and changes there probable need to be propagated here.

int NamdCentLB::requiredProxies(PatchID id, int neighborNodes[])
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

