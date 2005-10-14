
#if !defined(WIN32) || defined(__CYGWIN__)
#include <unistd.h>
#endif
#include <fcntl.h>

#include "InfoStream.h"
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

#if CHARM_VERSION > 050610
NamdCentLB::NamdCentLB(): CentralLB(CkLBOptions(-1))
#else
NamdCentLB::NamdCentLB()
#endif
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

  int numProcessors = count;
  int numPatches = PatchMap::Object()->numPatches();
  const int numComputes = ComputeMap::Object()->numComputes();
  const SimParameters* simParams = Node::Object()->simParameters;

  // these sizes should never change
  if ( ! processorArray ) processorArray = new processorInfo[numProcessors];
  if ( ! patchArray ) patchArray = new patchInfo[numPatches];
  if ( ! computeArray ) computeArray = new computeInfo[numComputes];

  int nMoveableComputes = buildData(stats,count);

  // gzheng debug
//#define DUMPDATA 1
//#define LOADDATA 1
#if DUMPDATA 
  dumpDataASCII("data", numProcessors, numPatches, nMoveableComputes);
#elif LOADDATA
  loadDataASCII("data", numProcessors, numPatches, nMoveableComputes);
//  dumpDataASCII("data.out", numProcessors, numPatches, nMoveableComputes);
//  CkExit();
#endif
  // end of debug section

  if (simParams->ldbStrategy == LDBSTRAT_REFINEONLY) {
    RefineOnly(computeArray,patchArray,processorArray,
                                nMoveableComputes, numPatches, numProcessors);
  } else if (simParams->ldbStrategy == LDBSTRAT_ALG7) {
    Alg7(computeArray,patchArray,processorArray,
                          nMoveableComputes, numPatches, numProcessors);
  } else if (simParams->ldbStrategy == LDBSTRAT_ALGORB) {
    if (step() == 1) {
      // iout << iINFO << "Load balance cycle " << step()
      //   << " using RecBisection\n" << endi;
      AlgRecBisection(computeArray,patchArray,processorArray,
                            nMoveableComputes, numPatches, numProcessors);
    } else {
      // iout << iINFO << "Load balance cycle " << step()
      //   << " using RefineOnly\n" << endi;
      RefineOnly(computeArray,patchArray,processorArray,
                                  nMoveableComputes, numPatches,
                                  numProcessors);
    }
  } else if (simParams->ldbStrategy == LDBSTRAT_OTHER) {
    // if (step() == 0) {
    if (step() < 2) {
      // iout << iINFO << "Load balance cycle " << step()
      //   << " using Alg7\n" << endi;
      Alg7(computeArray,patchArray,processorArray,
                            nMoveableComputes, numPatches, numProcessors);
    } else {
      // iout << iINFO << "Load balance cycle " << step()
      //   << " using RefineOnly\n" << endi;
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

#if LOADDATA
  dumpDataASCII("data.out", numProcessors, numPatches, nMoveableComputes);
  CkExit();
#endif

  // For error checking:
  // Count up computes, to see if somebody doesn't have any computes
  int i;
#if 0
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
#endif
  
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
#if CHARM_VERSION > 050611
  CLBMigrateMsg* msg = new(migrate_count,CkNumPes(),CkNumPes(),0) CLBMigrateMsg;
#else
  CLBMigrateMsg* msg = new(&migrate_count,1) CLBMigrateMsg;
#endif
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

void NamdCentLB::dumpDataASCII(char *file, int numProcessors,
			       int numPatches, int numComputes)
{
  char filename[128];
  sprintf(filename, "%s.%d", file, step());
  FILE* fp = fopen(filename,"w");
  if (fp == NULL){
     perror("dumpLDStatsASCII");
     return;
  }
  CkPrintf("***** DUMP data to file: %s ***** \n", filename);
  fprintf(fp,"%d %d %d\n",numProcessors,numPatches,numComputes);

  int i;
  for(i=0;i<numProcessors;i++) {
    processorInfo* p = processorArray + i;
    fprintf(fp,"%d %e %e %e %e\n",p->Id,p->load,p->backgroundLoad,p->computeLoad,p->idleTime);
  }

  for(i=0;i < numPatches; i++) {
    patchInfo* p = patchArray + i;
    fprintf(fp,"%d %e %d %d\n",p->Id,p->load,p->processor,p->numAtoms);
  }
    
  for(i=0; i < numComputes; i++) {
    computeInfo* c = computeArray + i;
    fprintf(fp,"%d %e %d %d %d %d",c->Id,c->load,c->patch1,c->patch2,
	    c->processor,c->oldProcessor);
#if CHARM_VERSION > 50910
    fprintf(fp, " %e %e", c->minTime, c->maxTime);
#endif
    fprintf(fp, "\n");
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
  CkExit();
}

void NamdCentLB::loadDataASCII(char *file, int &numProcessors,
			       int &numPatches, int &numComputes)
{
  char filename[128];
  sprintf(filename, "%s.%d", file, step());

  CkPrintf("***** Load ascii data from file: %s ***** \n", filename);

  FILE* fp = fopen(filename, "r");
  if (fp == NULL){
     perror("loadDataASCII");
     return;
  }

  fscanf(fp,"%d %d %d",&numProcessors,&numPatches,&numComputes);

  printf("numProcs: %d numPatches: %d numComputes: %d\n", numProcessors,numPatches, numComputes);

  delete [] processorArray;
  delete [] patchArray;
  delete [] computeArray;
  processorArray = new processorInfo[numProcessors];
  patchArray = new patchInfo[numPatches];
  computeArray = new computeInfo[numComputes];

  int i;
  for(i=0;i<numProcessors;i++) {
    processorInfo* p = processorArray + i;
    fscanf(fp,"%d %le %le %le", &p->Id, &p->load, &p->backgroundLoad, &p->computeLoad);
    fscanf(fp,"%le\n", &p->idleTime);
    if (p->Id != i) CmiAbort("Reading processorArray error!");
//    p->backgroundLoad = 0.0;
  }

  for(i=0;i < numPatches; i++) {
    patchInfo* p = patchArray + i;
    fscanf(fp,"%d %le %d %d\n",&p->Id,&p->load,&p->processor,&p->numAtoms);
    if (p->Id != i || p->processor > numProcessors || p->processor < 0) 
      CmiAbort("Reading patchArray error!");
  }
    
  for(i=0; i < numComputes; i++) {
    computeInfo* c = computeArray + i;
    fscanf(fp,"%d %le %d %d %d %d",&c->Id,&c->load,&c->patch1,&c->patch2,
	    &c->processor,&c->oldProcessor);
#if CHARM_VERSION > 50910
    fscanf(fp, " %le %le", &c->minTime, &c->maxTime);
#endif
    if (c->patch1 < 0 || c->patch1 > numPatches || c->patch2 < 0 || c->patch2 > numPatches)
      CmiAbort("Reading computeArray error!");
//  printf("%d %e %d %d %d %d %e %e\n", c->Id,c->load,c->patch1,c->patch2,c->processor,c->oldProcessor,c->minTime,c->maxTime);
  }

  // dump patchSet
  for (i=0; i< numProcessors; i++) {
      int num = processorArray[i].proxies.numElements();
      fscanf(fp,"%d",&num);
      for (int j=0; j<num; j++) {
          int id;
          fscanf(fp,"%d",&id);
          processorArray[i].proxies.insert(&patchArray[id]);
      }
  }
  // dump proxiesOn
  for (i=0; i<numPatches; i++)  {
      int num;
      fscanf(fp,"%d",&num);
      for (int j=0; j<num; j++) {
          int id;
	  fscanf(fp,"%d",&id);
          patchArray[i].proxiesOn.insert(&processorArray[id]);
      }
  }

  fclose(fp);
}
#endif

extern int isPmeProcessor(int); 

int NamdCentLB::buildData(CentralLB::LDStats* stats, int count)
{
  PatchMap* patchMap = PatchMap::Object();
  ComputeMap* computeMap = ComputeMap::Object();
  const SimParameters* simParams = Node::Object()->simParameters;

  BigReal bgfactor = simParams->ldbBackgroundScaling;
  BigReal pmebgfactor = simParams->ldbPMEBackgroundScaling;
  BigReal homebgfactor = simParams->ldbHomeBackgroundScaling;
  int pmeOn = simParams->PMEOn;
  int unLoadPme = simParams->ldbUnloadPME;
  int pmeBarrier = simParams->PMEBarrier;
  int unLoadZero = simParams->ldbUnloadZero;
  int unLoadRankZero = simParams->ldbUnloadRankZero;
  int unLoadSMP = simParams->ldbUnloadSMP;

  int i;
  for (i=0; i<count; ++i) {
    processorArray[i].Id = i;
    processorArray[i].available = CmiTrue;
    if ( pmeOn && isPmeProcessor(i) ) {
#if CHARM_VERSION > 050607
      processorArray[i].backgroundLoad = pmebgfactor * stats->procs[i].bg_walltime;
#else
      processorArray[i].backgroundLoad = pmebgfactor * stats[i].bg_walltime;
#endif
    } else if (patchMap->numPatchesOnNode(i) > 0) {
#if CHARM_VERSION > 050607
      processorArray[i].backgroundLoad = homebgfactor * stats->procs[i].bg_walltime;
#else
      processorArray[i].backgroundLoad = homebgfactor * stats[i].bg_walltime;
#endif
    } else {
#if CHARM_VERSION > 050607
      processorArray[i].backgroundLoad = bgfactor * stats->procs[i].bg_walltime;
#else
      processorArray[i].backgroundLoad = bgfactor * stats[i].bg_walltime;
#endif
    }
    processorArray[i].idleTime = stats->procs[i].idletime;
  }

#if 0
  double bgfactor = 1.0 + 1.0 * CkNumPes()/1000.0;
  if ( bgfactor > 2.0 ) bgfactor = 2.0;
  iout << iINFO << "Scaling background load by " << bgfactor << ".\n" << endi;
  int i;
  for (i=0; i<count; i++) {
    processorArray[i].Id = i;
    processorArray[i].backgroundLoad = bgfactor * stats[i].bg_walltime;
  }

  double bg_weight = 0.7;

  int i;
  for (i=0; i<count; i++) {
    processorArray[i].Id = i;
    if (patchMap->numPatchesOnNode(i) > 0)
#if CHARM_VERSION > 050607
      processorArray[i].backgroundLoad = bg_weight * stats->procs[i].bg_walltime;
#else
      processorArray[i].backgroundLoad = bg_weight * stats[i].bg_walltime;
#endif
    else 
#if CHARM_VERSION > 050607
      processorArray[i].backgroundLoad = stats[i].bg_walltime;
#else
      processorArray[i].backgroundLoad = stats->procs[i].bg_walltime;
#endif
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

  if (unLoadZero) processorArray[0].available = CmiFalse;
  if (unLoadRankZero) {
    for (int i=0; i<count; i+=4) 
      processorArray[i].available = CmiFalse;
  }

  // if all pes are Pme, disable this flag
  if (pmeOn && unLoadPme) {
    for (i=0; i<count; i++) {
      if (!isPmeProcessor(i))  break;
    }
    if (i==count) {
      iout << iINFO << "Turned off unLoadPme flag!\n"  << endi;
      unLoadPme = 0;
    }
  }
  
  if (pmeOn && unLoadPme) {
    for (i=0; i<count; i++) {
      if ((pmeBarrier && i==0) || isPmeProcessor(i)) 
	processorArray[i].available = CmiFalse;
    }
  }

  if (unLoadSMP) {
    int ppn = simParams->procsPerNode;
    int unloadrank = simParams->ldbUnloadRank;
    for (int i=0; i<count; i+=ppn) {
      processorArray[i+unloadrank].available = CmiFalse;
    }
  }

  int nMoveableComputes=0;
  int nProxies = 0;		// total number of estimated proxies
#if CHARM_VERSION > 050607
  int j;
  for (j=0; j < stats->n_objs; j++) {
      const LDObjData &this_obj = stats->objData[j];
      int frompe = stats->from_proc[j];
#else
  for (i=0; i < count; i++) {
    int j;
    for (j=0; j < stats[i].n_objs; j++) {
      const LDObjData &this_obj = stats[i].objData[j];
      int frompe = i;
#endif
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
#if CHARM_VERSION > 050607
	patchArray[pid].processor = stats->from_proc[j];
#else
	patchArray[pid].processor = i;
#endif
	const int numProxies = 
#if CMK_VERSION_BLUEGENE
	  requiredProxiesOnProcGrid(pid,neighborNodes);
#else
	requiredProxies(pid, neighborNodes);
#endif

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
#if CHARM_VERSION > 050607
	computeArray[nMoveableComputes].oldProcessor = stats->from_proc[j];
#else
	computeArray[nMoveableComputes].oldProcessor = i;
#endif
	computeArray[nMoveableComputes].processor = -1;
	computeArray[nMoveableComputes].patch1 = p0;
	computeArray[nMoveableComputes].patch2 = p1;
	computeArray[nMoveableComputes].handle = this_obj.handle;
	computeArray[nMoveableComputes].load = this_obj.wallTime;
#if CHARM_VERSION > 50910
	computeArray[nMoveableComputes].minTime = this_obj.minWall;
	computeArray[nMoveableComputes].maxTime = this_obj.maxWall;
#endif
	nMoveableComputes++;
      }
    }
#if ! ( CHARM_VERSION > 050607 )
  }
#endif
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
  stats->clear();
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
  PatchID neighbors[1 + PatchMap::MaxOneAway + PatchMap::MaxTwoAway];

  PatchMap* patchMap = PatchMap::Object();

  int myNode = patchMap->node(id);
  neighbors[0] = id;
  int numNeighbors = 1 + patchMap->downstreamNeighbors(id,neighbors+1);
  for ( i = 0; i < numNeighbors; ++i )
  {
    const int proxyNode = patchMap->basenode(neighbors[i]);
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
  int emptyNodes = numNodes - numPatches;
  if ( emptyNodes > numPatches ) {
    int nodesPerPatch = nProxyNodes + 1 + (emptyNodes-1) / numPatches;
    int proxyNode = (myNode + 1) % numNodes;
    while ( nProxyNodes < nodesPerPatch &&
			! patchMap->numPatchesOnNode(proxyNode) ) {
      if (proxyNode != myNode && proxyNodes[proxyNode] == No) {
        proxyNodes[proxyNode] = Yes;
        neighborNodes[nProxyNodes] = proxyNode;
        nProxyNodes++;
      }
      proxyNode = (proxyNode + 1) % numNodes;
    }
    proxyNode = (myNode - 1 + numNodes) % numNodes;
    while ( nProxyNodes < nodesPerPatch &&
			! patchMap->numPatchesOnNode(proxyNode) ) {
      if (proxyNode != myNode && proxyNodes[proxyNode] == No) {
        proxyNodes[proxyNode] = Yes;
        neighborNodes[nProxyNodes] = proxyNode;
        nProxyNodes++;
      }
      proxyNode = (proxyNode - 1 + numNodes) % numNodes;
    }
    proxyNode = (myNode + 1) % numNodes;
    while ( nProxyNodes < nodesPerPatch ) {
      if ( ! patchMap->numPatchesOnNode(proxyNode) &&
           proxyNode != myNode && proxyNodes[proxyNode] == No) {
        proxyNodes[proxyNode] = Yes;
        neighborNodes[nProxyNodes] = proxyNode;
        nProxyNodes++;
      }
      proxyNode = (proxyNode + 1) % numNodes;
    }
  } else {
    int proxyNode = myNode - 1;
    if ( proxyNode >= 0 && ! patchMap->numPatchesOnNode(proxyNode) ) {
      if (proxyNode != myNode && proxyNodes[proxyNode] == No) {
        proxyNodes[proxyNode] = Yes;
        neighborNodes[nProxyNodes] = proxyNode;
        nProxyNodes++;
      }
    }
    proxyNode = myNode + 1;
    if ( proxyNode < numNodes && ! patchMap->numPatchesOnNode(proxyNode) ) {
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

#if CMK_VERSION_BLUEGENE
// Figure out which proxies we will definitely create on other nodes,
// without regard for non-bonded computes.  This code is swiped from
// ProxyMgr, and changes there probable need to be propagated here.
// The proxies are placed on nearby processors on the 3d-grid along
// the X,Y,Z dimensions

 int NamdCentLB::requiredProxiesOnProcGrid(PatchID id, int neighborNodes[])
{
  enum proxyHere { No, Yes };
  int numNodes = CkNumPes();
  proxyHere *proxyNodes = new proxyHere[numNodes];
  int nProxyNodes;
  int i,j,k;

  int xsize = 0, ysize = 0, zsize = 0;
  int my_x =0, my_y = 0, my_z = 0;

  PatchMap* patchMap = PatchMap::Object();
  int myNode = patchMap->node(id);
    
  BGLTorusManager *tmanager = BGLTorusManager::getObject();
  xsize = tmanager->getXSize();
  ysize = tmanager->getYSize();
  zsize = tmanager->getZSize();
  
  tmanager->getCoordinatesByRank(myNode, my_x, my_y, my_z);
  
  if(xsize * ysize * zsize != CkNumPes()) {
    delete [] proxyNodes;
    return requiredProxies(id, neighborNodes);
  }  


  // Note all home patches.
  for ( i = 0; i < numNodes; ++i )
  {
    proxyNodes[i] = No;
  }
  nProxyNodes=0;

  // Check all two-away neighbors.
  // This is really just one-away neighbors, since 
  // two-away always returns zero: RKB
  PatchID neighbors[1 + PatchMap::MaxOneAway + PatchMap::MaxTwoAway];

  //Assign a proxy to all your neighbors. But dont increment counter
  //because these have to be there anyway.
  
  neighbors[0] = id;
  int numNeighbors = 1 + patchMap->downstreamNeighbors(id,neighbors+1);
  for ( i = 0; i < numNeighbors; ++i )
  {
    int proxyNode = patchMap->basenode(neighbors[i]);
    
    if (proxyNode != myNode)
      if (proxyNodes[proxyNode] == No)
	{
	  proxyNodes[proxyNode] = Yes;
	  neighborNodes[nProxyNodes] = proxyNode;
	  nProxyNodes++;
	}
    /*
    int px, py, pz;
    tmanager->getCoordinatesByRank(proxyNode, px, py, pz);
    
    //Place proxy in the mid point processor
    proxyNode = tmanager->coords2rank((my_x+px)/2, (my_y+py)/2, (my_z+pz)/2);
    
    if (proxyNode != myNode)
      if (proxyNodes[proxyNode] == No)
	{
	  proxyNodes[proxyNode] = Yes;
	  neighborNodes[nProxyNodes] = proxyNode;
	  nProxyNodes++;
	}
    */
  }
  
  //Place numNodesPerPatch proxies on the 3d torus neighbors of a processor

  int numPatches = patchMap->numPatches();
  int emptyNodes = numNodes - numPatches;
  //if ( emptyNodes > numPatches ) {
  
  int nodesPerPatch = nProxyNodes + 8 * (emptyNodes-1) / numPatches + 1;
  int proxyNode = 0 ;
  int proxy_x=0, proxy_y=0, proxy_z=0;
  
  //Choose from the 26 neighbors of mynode.
  //CkAssert(nodesPerPatch - nProxyNodes <= 26);  
  //Too few patches otherwise, try twoaway?
  
  for(k=-1; k<= 1; k++) {
    proxy_z = (my_z + k + zsize) % zsize;
    for(j=-1; j <= 1; j++) {
      proxy_y = (my_y + j + ysize) % ysize;
      for(i = -1; i <= 1; i++) {
	if(i == 0 && j == 0 && k == 0)
	  continue;
	
	proxy_x = (my_x + i + xsize) % xsize;
	proxyNode = tmanager->coords2rank(proxy_x, proxy_y, proxy_z);

	if(! patchMap->numPatchesOnNode(proxyNode) &&
	   proxyNodes[proxyNode] == No) {
	  proxyNodes[proxyNode] = Yes;
	  neighborNodes[nProxyNodes] = proxyNode;
	  nProxyNodes++;
	}
	
	if(nProxyNodes >= nodesPerPatch || 
	   nProxyNodes >= PatchMap::MaxOneAway + PatchMap::MaxTwoAway)
	  break;	  
      }
      
      if(nProxyNodes >= nodesPerPatch || 
	 nProxyNodes >= PatchMap::MaxOneAway + PatchMap::MaxTwoAway)
	break;	  
    }
    if(nProxyNodes >= nodesPerPatch || 
       nProxyNodes >= PatchMap::MaxOneAway + PatchMap::MaxTwoAway)
      break;	  
  }        
  //  } 
  //else
  //CkAbort("NumPes < 2*numpatches\n\n");

  CkPrintf("Returning %d proxies\n", nProxyNodes);

  delete [] proxyNodes;
  return nProxyNodes;
}

#endif
