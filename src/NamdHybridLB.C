/*****************************************************************************
 * $Source: /home/cvs/namd/cvsroot/namd2/src/NamdHybridLB.C,v $
 * $Author: jim $
 * $Date: 2010/03/08 17:34:21 $
 * $Revision: 1.5 $
 *****************************************************************************/

#if !defined(WIN32) || defined(__CYGWIN__)
#include <unistd.h>
#endif
#include <fcntl.h>

#include "InfoStream.h"
#include "NamdHybridLB.h"
#include "NamdHybridLB.def.h"
#include "Node.h"
#include "PatchMap.h"
#include "ComputeMap.h"
#include "LdbCoordinator.h"

extern int isPmeProcessor(int); 

/**
 * Creates the chare array for the hybrid load balancer.
 */ 
void CreateNamdHybridLB() {
  CProxy_NamdHybridLB::ckNew();
}

/**
 * @brief Default constructor.
 */
NamdHybridLB::NamdHybridLB(): HybridBaseLB(CkLBOptions(-1))
{

	// setting the name
	lbname = (char *)"NamdHybridLB";

	// initializing thisProxy
	thisProxy = CProxy_NamdHybridLB(thisgroup);
	
	// initializing the central LB
	centralLB = AllocateNamdCentLB();

	// initializing the dummy LB
	dummyLB = AllocateNamdDummyLB();

	// assigning initial values to variables
	computeArray = NULL;
	patchArray = NULL;
	processorArray = NULL;
	updateCount = 0;
	updateFlag = false;
	collectFlag = false;

}

/**
 * @brief Function used to discover if load balancer can be run at that point.
 *
 * It is called from HybridBase every time AtSync method is called.
 */
CmiBool NamdHybridLB::QueryBalanceNow(int _step){ 
	if ( LdbCoordinator::Object()->takingLdbData ) {
		return CmiTrue;
	} else {
		return CmiFalse;
	} 
}

CmiBool NamdHybridLB::QueryDumpData(){                                                                                                 
#if 0                                                                                             
  if (LdbCoordinator::Object()->ldbCycleNum == 1)  return CmiTrue;                                
  if (LdbCoordinator::Object()->ldbCycleNum == 2)  return CmiTrue;                                
#endif                                                                                            
  return CmiFalse;                                                                                
}

/**
 *  * Runs the load balancing strategy with shrinking the load information.
 *   * Note: now, it is just calling the Strategy, but eventually will have
 *    * its own code.
 *     */
LBVectorMigrateMsg* NamdHybridLB::VectorStrategy(LDStats* stats,int count){
		CkPrintf("[%d] Using Vector Strategy to balance the load\n",CkMyPe());
        LBVectorMigrateMsg* msg = new(0,0) LBVectorMigrateMsg;
        msg->n_moves = 0;
        msg->level = currentLevel;
        return msg;
}

/*
 * Runs the load balancing strategy
 */
CLBMigrateMsg* NamdHybridLB::Strategy(LDStats* stats,int count){
	
  	// CkPrintf("[%d] NamdHybridLB at Strategy\n",CkMyPe());
	
	// calling the centralLB for level 1		
	if(currentLevel == 1){
		LevelData *lData = levelData[currentLevel];
		CLBMigrateMsg *msg, *newMsg;
		msg = CentralStrategy(stats,count);

		// creating a new message to send to its parent
		newMsg = new(msg->n_moves,CkNumPes(),CkNumPes(),0) CLBMigrateMsg;
  		newMsg->level = currentLevel;
  		newMsg->n_moves = msg->n_moves;
  		for(int i=0; i < msg->n_moves; i++) {
    		newMsg->moves[i] = msg->moves[i];
  		}
		thisProxy[0].UpdateComputeMap(newMsg);
		return msg;
	}else{
		dummyLB->work(stats,count);	
		return createMigrateMsg(stats,count);
	}
}


/**
 * Updates the compute map with the migration information from its children
 */
void NamdHybridLB::UpdateComputeMap(CLBMigrateMsg *msg){
	int children;

	// getting the number of children
	children = tree->numNodes(currentLevel);
	// CkPrintf("[%d] Updating compute map, total %d\n",CkMyPe(),siblings);

	// getting the compute map to insert the changes coming from the children
	ComputeMap *computeMap = ComputeMap::Object();

	// traversing the set of moves in msg
	for(int i=0; i<msg->n_moves; i++){
		computeMap->setNewNode(msg->moves[i].obj.id.id[0],msg->moves[i].to_pe);	
	}

	// checking if all children have sent the update
	updateCount++;
	if(updateCount == children){
		updateCount = 0;
		updateFlag = true;
		 // CkPrintf("[%d] UPDATE READY\n",CkMyPe());		
	}

	// checking if the collect info is ready
	if(updateFlag && collectFlag){
		updateFlag = false;
		collectFlag = false;	
		thisProxy[parent_backup].CollectInfo(loc_backup, n_backup, fromlevel_backup);
	}

}

/**
 * This function implements a strategy similar to the one used in the 
 * centralized case in NamdCentLB.
 */
CLBMigrateMsg* NamdHybridLB::CentralStrategy(LDStats* stats,int count){

  //  CkPrintf("LDB: All statistics received at %f, %f\n",
  //  CmiTimer(),CmiWallTimer());

  int numProcessors = count;
  int numPatches = PatchMap::Object()->numPatches();
  ComputeMap *computeMap = ComputeMap::Object();
  const int numComputes = computeMap->numComputes();
  const SimParameters* simParams = Node::Object()->simParameters;

  // these sizes should never change
  if ( ! processorArray ) processorArray = new processorInfo[numProcessors];
  if ( ! patchArray ) patchArray = new patchInfo[numPatches];
  if ( ! computeArray ) computeArray = new computeInfo[numComputes];


  int nMoveableComputes = buildData(stats,count);

#if LDB_DEBUG
#define DUMP_LDBDATA 1
#define LOAD_LDBDATA 1
#endif

#if DUMP_LDBDATA 
  dumpDataASCII("ldbd_before", numProcessors, numPatches, nMoveableComputes);
#elif LOAD_LDBDATA
  loadDataASCII("ldbd_before.5", numProcessors, numPatches, nMoveableComputes);
  // CkExit();
#endif

	if (step() < 2)
    	TorusLB(computeArray, patchArray, processorArray,
	          nMoveableComputes, numPatches, numProcessors);
    else
      	RefineTorusLB(computeArray, patchArray, processorArray,
                  nMoveableComputes, numPatches, numProcessors, 1);

#if LDB_DEBUG && USE_TOPOMAP
  TopoManager tmgr;
  int pe1, pe2, pe3, hops=0;
  /* This is double counting the hops
  for(int i=0; i<nMoveableComputes; i++)
  {
    pe1 = computeArray[i].processor;
    pe2 = patchArray[computeArray[i].patch1].processor;
    pe3 = patchArray[computeArray[i].patch2].processor;
    hops += tmgr.getHopsBetweenRanks(pe1, pe2);
    if(computeArray[i].patch1 != computeArray[i].patch2)
      hops += tmgr.getHopsBetweenRanks(pe1, pe3);  
  }*/
  for (int i=0; i<numPatches; i++)  {
    //int num = patchArray[i].proxiesOn.numElements();
    pe1 = patchArray[i].processor;
    Iterator nextProc;
    processorInfo *p = (processorInfo *)patchArray[i].proxiesOn.iterator((Iterator *)&nextProc);
    while (p) {
      pe2 = p->Id;
      hops += tmgr.getHopsBetweenRanks(pe1, pe2);
      p = (processorInfo *)patchArray[i].proxiesOn.next((Iterator*)&nextProc);
    }
  }
  CkPrintf("Load Balancing: Number of Hops: %d\n", hops);
#endif

#if DUMP_LDBDATA
  dumpDataASCII("ldbd_after", numProcessors, numPatches, nMoveableComputes);
#elif LOAD_LDBDATA
  dumpDataASCII("ldbd_after.5", numProcessors, numPatches, nMoveableComputes);
  // loadDataASCII("ldbd_after", numProcessors, numPatches, nMoveableComputes);
  // CkExit();
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
            /* ERASECkPrintf("[%d] Obj %d migrating from %d to %d\n",
                     CkMyPe(),computeArray[i].handle.id.id[0],
			 computeArray[i].oldProcessor, computeArray[i].processor); */
      MigrateInfo *migrateMe = new MigrateInfo;
      migrateMe->obj = computeArray[i].handle;
      migrateMe->from_pe = computeArray[i].oldProcessor;
      migrateMe->to_pe = computeArray[i].processor;
      migrateInfo.insertAtEnd(migrateMe);

      // sneak in updates to ComputeMap
      //ERASE CkPrintf("%d setting %d to processor %d\n",CkMyPe(),computeArray[i].handle.id.id[0],computeArray[i].processor);
      computeMap->setNewNode(computeArray[i].handle.id.id[0],
				computeArray[i].processor);
    }
  }
  // CkPrintf("LOAD BALANCING READY %d\n",CkMyPe()); 

 
  int migrate_count=migrateInfo.length();
  // CkPrintf("NamdCentLB migrating %d elements\n",migrate_count);
  CLBMigrateMsg* msg = new(migrate_count,CkNumPes(),CkNumPes(),0) CLBMigrateMsg;
  msg->level = currentLevel;
  msg->n_moves = migrate_count;
  for(i=0; i < migrate_count; i++) {
    MigrateInfo* item = migrateInfo[i];
    msg->moves[i] = *item;
    delete item;
    migrateInfo[i] = 0;
  }

/*Not needed  for (i=0; i<numProcessors; i++) {
    cpuloads[i] = processorArray[i].load;
  }
*/

  delete [] processorArray;
  delete [] patchArray;
  delete [] computeArray;

  processorArray = NULL;
  patchArray = NULL;
  computeArray = NULL;
  
  return msg;

}

/**
 * @brief Builds the data structures required for the load balancing strategies in NAMD.
 */ 
int NamdHybridLB::buildData(CentralLB::LDStats* stats, int count){
	int i;

   // CkPrintf("%d is in buildData %d %d %d\n",CmiMyPe(),stats->count,stats->complete_flag,stats->procs[1].pe);

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
	int unLoadOne = simParams->ldbUnloadOne;

	// traversing the list of processors and getting their load information
	for (i=0; i<count; ++i) {
		//BACKUPprocessorArray[i].Id = i; 
		processorArray[i].Id = stats->procs[i].pe;
    	processorArray[i].available = CmiTrue;
    	//BACKUPif ( pmeOn && isPmeProcessor(i) ) {
    	if ( pmeOn && isPmeProcessor(stats->procs[i].pe) ) {
     		processorArray[i].backgroundLoad = pmebgfactor * stats->procs[i].bg_walltime;
    	//BACKUP } else if (patchMap->numPatchesOnNode(i) > 0) {
    	} else if (patchMap->numPatchesOnNode(stats->procs[i].pe) > 0) {
      		processorArray[i].backgroundLoad = homebgfactor * stats->procs[i].bg_walltime;
    	} else {
      		processorArray[i].backgroundLoad = bgfactor * stats->procs[i].bg_walltime;
    	}
    	processorArray[i].idleTime = stats->procs[i].idletime;
    	processorArray[i].load = processorArray[i].computeLoad = 0.0;
	}


	if(unLoadZero && stats->procs[0].pe == 0) processorArray[0].available = CmiFalse;
	if(unLoadOne && stats->procs[1].pe == 1) processorArray[1].available = CmiFalse;

	// if all pes are Pme, disable this flag
  	if (pmeOn && unLoadPme) {
    	for (i=0; i<count; i++) {
      		if(!isPmeProcessor(stats->procs[i].pe))  break;
    	}
    	if (i==count) {
      		iout << iINFO << "Turned off unLoadPme flag!\n"  << endi;
      		unLoadPme = 0;
		}
  	}
  
	if (pmeOn && unLoadPme) {
    	for (i=0; i<count; i++) {
      		if ((pmeBarrier && i==0) || isPmeProcessor(stats->procs[i].pe)) 
				processorArray[i].available = CmiFalse;
    	}
  	}

	int nMoveableComputes=0;
	int nProxies = 0;		// total number of estimated proxies
	int index;

	int j;
  	for(j=0; j < stats->n_objs; j++){
		const LDObjData &this_obj = stats->objData[j];
      	int frompe = stats->from_proc[j];

		// filter out non-NAMD managed objects (like PME array)
      	if (this_obj.omID().id.idx != 1) continue;

      	if (this_obj.id().id[1] == -2) { // Its a patch
			const int pid = this_obj.id().id[0];
			int neighborNodes[PatchMap::MaxOneAway + PatchMap::MaxTwoAway];

			patchArray[pid].Id = pid;
			patchArray[pid].numAtoms = 0;
			patchArray[pid].processor = stats->from_proc[j];
			const int numProxies = 
#if 0 // USE_TOPOMAP - this function needs to be there for the hybrid case
			requiredProxiesOnProcGrid(pid,neighborNodes);
#else
			requiredProxies(pid, neighborNodes);
#endif

        	nProxies += numProxies;

			for (int k=0; k<numProxies; k++) {
				if( (neighborNodes[k] >= stats->procs[0].pe) && (neighborNodes[k] <= stats->procs[count-1].pe) ){
					index = neighborNodes[k] - stats->procs[0].pe;
	  				//BACKUP processorArray[neighborNodes[k]].proxies.insert(&patchArray[pid]);
	  				processorArray[index].proxies.insert(&patchArray[pid]);
	  				//BACKUP patchArray[pid].proxiesOn.insert(&processorArray[neighborNodes[k]]);
	  				patchArray[pid].proxiesOn.insert(&processorArray[index]);
				}
			}
      	} else if (this_obj.migratable) { // Its a compute

			const int cid = this_obj.id().id[0];
			const int p0 = computeMap->pid(cid,0);

			// For self-interactions, just return the same pid twice
			int p1;
			if (computeMap->numPids(cid) > 1)
	  			p1 = computeMap->pid(cid,1);
			else p1 = p0;
			computeArray[nMoveableComputes].Id = cid;
			//BACKUP computeArray[nMoveableComputes].oldProcessor = stats->from_proc[j];
			computeArray[nMoveableComputes].oldProcessor = stats->from_proc[j] + stats->procs[0].pe;

			index = stats->from_proc[j]; 
			//BACKUP2 index = stats->from_proc[j] - stats->procs[0].pe;
			//BACKUP processorArray[stats->from_proc[j]].computeLoad += this_obj.wallTime;
			processorArray[index].computeLoad += this_obj.wallTime;
			computeArray[nMoveableComputes].processor = -1;
			computeArray[nMoveableComputes].patch1 = p0;
			computeArray[nMoveableComputes].patch2 = p1;
			computeArray[nMoveableComputes].handle = this_obj.handle;
			computeArray[nMoveableComputes].load = this_obj.wallTime;
			computeArray[nMoveableComputes].minTime = this_obj.minWall;
			computeArray[nMoveableComputes].maxTime = this_obj.maxWall;
			nMoveableComputes++;
      	}
	}

  	for (i=0; i<count; i++) {
    	processorArray[i].load = processorArray[i].backgroundLoad + processorArray[i].computeLoad;
  	}
  	stats->clear();
  	return nMoveableComputes;
}


int NamdHybridLB::requiredProxies(PatchID id, int neighborNodes[])
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
  
  // This code needs to be turned off when the creation of ST is
  // shifted to the load balancers -ASB

#if 1
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
    int count = 0;
    while ( nProxyNodes < nodesPerPatch ) {
      if ( ! patchMap->numPatchesOnNode(proxyNode) &&
           proxyNode != myNode && proxyNodes[proxyNode] == No) {
        proxyNodes[proxyNode] = Yes;
        neighborNodes[nProxyNodes] = proxyNode;
        nProxyNodes++;
      }
      proxyNode = (proxyNode + 1) % numNodes;
      count ++; if (count == numNodes) break;   // we looped all
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
#endif

  delete [] proxyNodes;
  return nProxyNodes;
}

void NamdHybridLB::CollectInfo(Location *loc, int n, int fromlevel)
{
   int atlevel = fromlevel + 1;
   LevelData *lData = levelData[atlevel];
   lData->info_recved++;

   CkVec<Location> &matchedObjs = lData->matchedObjs;

   // sort into mactched and unmatched list
#if CHARM_VERSION <= 60103
   CkVec<Location> &unmatchedObjs = lData->unmatchedObjs;
   for (int i=0; i<n; i++) {
     // search and see if we have answer, put to matched
     // store in unknown
     int found = 0;
     for (int obj=0; obj<unmatchedObjs.size(); obj++) {
       if (loc[i].key == unmatchedObjs[obj].key) {
         // answer must exist
         CmiAssert(unmatchedObjs[obj].loc != -1 || loc[i].loc != -1);
         if (unmatchedObjs[obj].loc == -1) unmatchedObjs[obj].loc = loc[i].loc;
         matchedObjs.push_back(unmatchedObjs[obj]);
         unmatchedObjs.remove(obj);
         found = 1;
         break;
       }
     }
     if (!found) unmatchedObjs.push_back(loc[i]);
   }
#else
   std::map<LDObjKey, int> &unmatchedObjs = lData->unmatchedObjs;
   for (int i=0; i<n; i++) {
     std::map<LDObjKey, int>::iterator iter = unmatchedObjs.find(loc[i].key);
     if (iter != unmatchedObjs.end()) {
       CmiAssert(iter->second != -1 || loc[i].loc != -1);
       if (loc[i].loc == -1) loc[i].loc = iter->second;
       matchedObjs.push_back(loc[i]);
       unmatchedObjs.erase(iter);
     }
     else
       unmatchedObjs[loc[i].key] = loc[i].loc;
   }
#endif

//  DEBUGF(("[%d] level %d has %d unmatched and %d matched. \n", CkMyPe(), atlevel, unmatchedObjs.size(), matchedObjs.size()));

   if (lData->info_recved == lData->nChildren) {
     lData->info_recved = 0;
     if (_lb_args.debug() > 1)
         CkPrintf("[%d] CollectInfo at level %d started at %f\n",
	        CkMyPe(), atlevel, CkWallTimer());
     if (lData->parent != -1) {

		// NAMD specific
#if CHARM_VERSION <= 60103
#define unmatchedbuf unmatchedObjs
#else
		CkVec<Location> unmatchedbuf;
   		for(std::map<LDObjKey, int>::const_iterator it = unmatchedObjs.begin(); it != unmatchedObjs.end(); ++it){
    		unmatchedbuf.push_back(Location(it->first, it->second));
   		}
#endif
		// checking if update of ComputeMap is ready before calling parent
		if(CkMyPe() == 0){
			if(updateFlag){
				updateFlag = false;
				collectFlag = false;
				thisProxy[lData->parent].CollectInfo(unmatchedbuf.getVec(), unmatchedbuf.size(), atlevel);
			}else{
				CkPrintf("[%d] COMPUTEMAP UPDATE NOT READY\n",CkMyPe());
				collectFlag = true;
				parent_backup = lData->parent;
				loc_backup = unmatchedbuf.getVec();
				n_backup = unmatchedbuf.size();
				fromlevel_backup = atlevel;
			}
		}else{
			// send only unmatched ones up the tree
			thisProxy[lData->parent].CollectInfo(unmatchedbuf.getVec(), unmatchedbuf.size(), atlevel);
		}

     }
     else { // root
       // we should have all answers now
       CmiAssert(unmatchedObjs.size() == 0);
       // start send match list down
       thisProxy.PropagateInfo(matchedObjs.getVec(), matchedObjs.size(), atlevel, lData->nChildren, lData->children);
       lData->statsData->clear();
     }
   }
}
