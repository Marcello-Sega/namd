
#if !defined(WIN32) || defined(__CYGWIN__)
#include <unistd.h>
#endif
#include <fcntl.h>

#include <charm++.h>

#include "NamdCentLB.h"
#include "NamdCentLB.def.h"
#include "Node.h"
#include "PatchMap.h"
#include "ComputeMap.h"
#include "LdbCoordinator.h"

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
  processorArray = 0;
  patchArray = 0;
  computeArray = 0;
}

/*
NamdCentLB::~NamdCentLB()
{
  delete [] processorArray;
  delete [] patchArray;
  delete [] computeArray;
}
*/

CmiBool NamdCentLB::QueryBalanceNow(int _step)
{
  //  CkPrintf("[%d] Balancing on step %d\n",CkMyPe(),_step);
  if ( LdbCoordinator::Object()->takingLdbData ) {
    return CmiTrue;
  } else {
    return CmiFalse;
  }
}

CmiBool NamdCentLB::QueryDumpData()
{
#if 0
  if (LdbCoordinator::Object()->ldbCycleNum == 1)  return CmiTrue;
  if (LdbCoordinator::Object()->ldbCycleNum == 2)  return CmiTrue;
#endif
  return CmiFalse;
}
            
CLBMigrateMsg* NamdCentLB::Strategy(CentralLB::LDStats* stats, int count)
{
  //  CkPrintf("LDB: All statistics received at %f, %f\n",
  //  CmiTimer(),CmiWallTimer());

  const int numProcessors = count;
  const int numPatches = PatchMap::Object()->numPatches();
  const int numComputes = ComputeMap::Object()->numComputes();
  const SimParameters* simParams = Node::Object()->simParameters;

  // these sizes should never change
  if ( ! processorArray ) processorArray = new processorInfo[numProcessors];
  if ( ! patchArray ) patchArray = new patchInfo[numPatches];
  if ( ! computeArray ) computeArray = new computeInfo[numComputes];

  const int nMoveableComputes = buildData(stats,count);
  // gzheng debug
  // dumpData("data", numProcessors, numPatches, nMoveableComputes);
  // loadData("data", numProcessors, numPatches, nMoveableComputes);
  // end of debug section

  if (simParams->ldbStrategy == LDBSTRAT_REFINEONLY) {
    RefineOnly(computeArray,patchArray,processorArray,
                                nMoveableComputes, numPatches, numProcessors);
  } else if (simParams->ldbStrategy == LDBSTRAT_ALG7) {
    Alg7(computeArray,patchArray,processorArray,
                          nMoveableComputes, numPatches, numProcessors);
  } else if (simParams->ldbStrategy == LDBSTRAT_ALGORB) {
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

  // For error checking:
  // Count up computes, to see if somebody doesn't have any computes
  int i;
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
  
  CkVec<MigrateInfo *> migrateInfo;
  for(i=0;i<nMoveableComputes;i++) {
    if (computeArray[i].processor != computeArray[i].oldProcessor) {
      //      CkPrintf("[%d] Obj %d migrating from %d to %d\n",
      //               CkMyPe(),computeArray[i].handle.id.id[0],
      //	       computeArray[i].processor,computeArray[i].oldProcessor);
      MigrateInfo *migrateMe = new MigrateInfo;
      migrateMe->obj = computeArray[i].handle;
      migrateMe->from_pe = computeArray[i].oldProcessor;
      migrateMe->to_pe = computeArray[i].processor;
      migrateInfo.insertAtEnd(migrateMe);
    }
  }
  
  int migrate_count=migrateInfo.length();
  // CkPrintf("NamdCentLB migrating %d elements\n",migrate_count);
  CLBMigrateMsg* msg = new(&migrate_count,1) CLBMigrateMsg;
  msg->n_moves = migrate_count;
  for(i=0; i < migrate_count; i++) {
    MigrateInfo* item = migrateInfo[i];
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
      int num = processorArray[i].proxies.numElements();
      write(fd, &num, sizeof(int));
//CkPrintf("**** Proc:%d num:%d \n", i, num);
      Iterator nextProxy;
//      nextProxy.id = 0;
      patchInfo *p = (patchInfo *)processorArray[i].proxies.iterator((Iterator *)&nextProxy);
      while (p) 
      {
          write(fd, &p->Id, sizeof(int));
          p = (patchInfo *)processorArray[i].proxies.next((Iterator*)&nextProxy);
      }
  }
  // dump proxiesOn
  for (i=0; i<numPatches; i++)
  {
      int num = patchArray[i].proxiesOn.numElements();
      write(fd, &num, sizeof(int));
// CkPrintf("**** Patch:%d num:%d \n", i, num);
      Iterator nextProc;
//      nextProc.id = 0;
      processorInfo *p = (processorInfo *)patchArray[i].proxiesOn.iterator((Iterator *)&nextProc);
      while (p) 
      {
          write(fd, &p->Id, sizeof(int));
          p = (processorInfo *)patchArray[i].proxiesOn.next((Iterator*)&nextProc);
      }
  }

  close(fd);
}

void NamdCentLB::dumpDataASCII(char *file, int numProcessors,
			       int numPatches, int numComputes)
{
  FILE* fp = fopen(file,"w");
  if (fp == NULL){
     perror("dumpLDStats");
     return;
  }
  CkPrintf("***** DUMP data to file: %s ***** \n", file);
  fprintf(fp,"%d %d %d\n",numProcessors,numPatches,numComputes);

  int i;
  for(i=0;i<numProcessors;i++) {
    processorInfo* p = processorArray + i;
    fprintf(fp,"%d %e %e %e\n",p->Id,p->load,p->backgroundLoad,p->computeLoad);
  }

  for(i=0;i < numPatches; i++) {
    patchInfo* p = patchArray + i;
    fprintf(fp,"%d %e %d %d\n",p->Id,p->load,p->processor,p->numAtoms);
  }
    
  for(i=0; i < numComputes; i++) {
    computeInfo* c = computeArray + i;
    fprintf(fp,"%d %e %d %d %d %d\n",c->Id,c->load,c->patch1,c->patch2,
	    c->processor,c->oldProcessor);
  }

  // dump patchSet
  for (i=0; i< numProcessors; i++) {
      int num = processorArray[i].proxies.numElements();
      fprintf(fp,"%d\n",num);
      Iterator nextProxy;
      patchInfo *p = (patchInfo *)processorArray[i].proxies.
	iterator((Iterator *)&nextProxy);
      while (p) {
          fprintf(fp,"%d\n",p->Id);
          p = (patchInfo *)processorArray[i].proxies.
	    next((Iterator*)&nextProxy);
      }
  }
  // dump proxiesOn
  for (i=0; i<numPatches; i++)  {
    int num = patchArray[i].proxiesOn.numElements();
    fprintf(fp,"%d\n",num);
      Iterator nextProc;
      processorInfo *p = (processorInfo *)patchArray[i].proxiesOn.
	iterator((Iterator *)&nextProc);
      while (p) {
	fprintf(fp,"%d\n",p->Id);
	p = (processorInfo *)patchArray[i].proxiesOn.
	  next((Iterator*)&nextProc);
      }
  }

  fclose(fp);
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

/*
  processorInfo *processors = new processorInfo[numProcessors];
  computeInfo *computes = new computeInfo[numComputes];
  patchInfo *patchs = new patchInfo [numPatches];
*/

  read(fd, processorArray, sizeof(processorInfo)*numProcessors);
  read(fd, patchArray, sizeof(patchInfo)*numPatches);
  read(fd, computeArray, sizeof(computeInfo)*numComputes);

  int i;
  for (i=0; i<numProcessors; i++)
  {
      int num;
      read(fd, &num, sizeof(int));
      printf("proc: %d proxies:%d \n", i, num);
      for (int j=0; j<num; j++) {
          int id;
          read(fd, &id, sizeof(int));
	  printf("%d ", id);
          processorArray[i].proxies.insert(&patchArray[id]);
      }
      printf("\n");
  }
  for (i=0; i<numPatches; i++)
  {
      int num;
      read(fd, &num, sizeof(int));
//      printf("patch: %d proxiesOn:%d \n", i, num);
      for (int j=0; j<num; j++) {
          int id;
          read(fd, &id, sizeof(int));
          patchArray[i].proxiesOn.insert(&processorArray[id]);
      }
  }

  close(fd);
}

#endif

extern int isPmeProcessor(int); 

int NamdCentLB::buildData(CentralLB::LDStats* stats, int count)
{
  PatchMap* patchMap = PatchMap::Object();
  ComputeMap* computeMap = ComputeMap::Object();

#if 1
  double bgfactor = 1.0 + 1.0 * CkNumPes()/1000.0;
  if ( bgfactor > 2.0 ) bgfactor = 2.0;
  iout << iINFO << "Scaling background load by " << bgfactor << ".\n" << endi;
  int i;
  for (i=0; i<count; i++) {
    processorArray[i].Id = i;
    processorArray[i].backgroundLoad = bgfactor * stats[i].bg_walltime;
  }

#else
  double bg_weight = 0.7;

  int i;
  for (i=0; i<count; i++) {
    processorArray[i].Id = i;
    if (patchMap->numPatchesOnNode(i) > 0)
      processorArray[i].backgroundLoad = bg_weight * stats[i].bg_walltime;
    else 
      processorArray[i].backgroundLoad = stats[i].bg_walltime;
  }
  
  //Modification to reduce the coputeload on PME processors
  const SimParameters* simParams = Node::Object()->simParameters;  
  
  // CkPrintf("BACKGROUND LOAD\n");
  if(simParams->PMEOn) {
    double bgfactor = 1.0 + 1.0 * CkNumPes()/1000.0;
    if ( bgfactor > 2.0 ) bgfactor = 2.0;
    for (i=0; i<count; i++) {
      // CkPrintf("BG[%d] =  %5.5lf,", i, processorArray[i].backgroundLoad);
      if(isPmeProcessor(i)) {
	processorArray[i].backgroundLoad *= bgfactor;
      }
      // CkPrintf("%5.5lf;  ", processorArray[i].backgroundLoad);
    }
  }
  // CkPrintf("\n");
#endif  

  int nMoveableComputes=0;
  int nProxies = 0;		// total number of estimated proxies
  for (i=0; i < count; i++) {
    int j;
    for (j=0; j < stats[i].n_objs; j++) {
      const LDObjData &this_obj = stats[i].objData[j];
      // filter out non-NAMD managed objects (like PME array)
#if CHARM_VERSION > 050405
      if (this_obj.omID().id.idx != 1) continue;
#elif CHARM_VERSION > 050403
      if (this_obj.omID.id.idx != 1) continue;
#else
      if (this_obj.omID.id != 1) continue;
#endif
#if CHARM_VERSION > 050405
      if (this_obj.id().id[1] == -2) { // Its a patch
	const int pid = this_obj.id().id[0];
#else
      if (this_obj.id.id[1] == -2) { // Its a patch
	const int pid = this_obj.id.id[0];
#endif
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
      } else if (this_obj.migratable) { // Its a compute
#if CHARM_VERSION > 050405
	const int cid = this_obj.id().id[0];
#else
	const int cid = this_obj.id.id[0];
#endif
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
#if 0
  int averageProxy = nProxies / count;
  CkPrintf("total proxies: %d, avervage: %d\n", nProxies, averageProxy);
  for (i=0; i<count; i++) {
    // too many proxies on this node, weight the background load
    int proxies = processorArray[i].proxies.numElements();
    if (proxies > averageProxy) {
      double factor = 1.0*(proxies-averageProxy)/nProxies;
      processorArray[i].backgroundLoad *= (1.0 + factor);
      CkPrintf("On [%d]: too many proxies: %d, increased bg load by %f\n", i, nProxies, factor);
    }
  }
#endif
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

  // Distribute initial default proxies across empty processors.
  // This shouldn't be necessary, but may constrain the load balancer
  // and avoid placing too many proxies on a single processor.  -JCP

  int numPatches = patchMap->numPatches();
  if ( numNodes > 1.25 * numPatches ) {  // avoid marginal cases
    int emptyNodes = numNodes - numPatches;
    int nodesPerPatch = 3 + emptyNodes / numPatches;
    int baseNode = (myNode - 1 + numNodes) % numNodes;
    for ( i = 0; i < nodesPerPatch; ++i ) {
      int proxyNode = (baseNode+i) % numNodes;
      if (proxyNode != myNode && proxyNodes[proxyNode] == No) {
	proxyNodes[proxyNode] = Yes;
	neighborNodes[nProxyNodes] = proxyNode;
	nProxyNodes++;
      }
    }
  }

  delete [] proxyNodes;
  return nProxyNodes;
}

