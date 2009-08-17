
/**
***  Copyright (c) 1995, 1996, 1997, 1998, 1999, 2000 by
***  The Board of Trustees of the University of Illinois.
***  All rights reserved.
**/

/*
   HomePatch owns the actual atoms of a Patch of space
   Proxy(s) get messages via ProxyMgr from HomePatch(es)
   to update lists of atoms and their coordinates
   HomePatch(es) also have a Sequencer bound to them

   superclass: 	Patch		
*/

#include <math.h>
#include "charm++.h"

#include "SimParameters.h"
#include "HomePatch.h"
#include "AtomMap.h"
#include "Node.h"
#include "PatchMap.inl"
#include "main.h"
#include "ProxyMgr.decl.h"
#include "ProxyMgr.h"
#include "Migration.h"
#include "Molecule.h"
#include "PatchMgr.h"
#include "Sequencer.h"
#include "LdbCoordinator.h"
#include "Settle.h"
#include "ReductionMgr.h"
#include "Sync.h"
#include "Random.h"
#include "Priorities.h"

#define TINY 1.0e-20;
#define MAXHGS 10
#define MIN_DEBUG_LEVEL 2
//#define DEBUGM
#include "Debug.h"

typedef int HGArrayInt[MAXHGS];
typedef BigReal HGArrayBigReal[MAXHGS];
typedef zVector HGArrayVector[MAXHGS];
typedef BigReal HGMatrixBigReal[MAXHGS][MAXHGS];
typedef zVector HGMatrixVector[MAXHGS][MAXHGS];

int average(CompAtom *qtilde,const HGArrayVector &q,BigReal *lambda,const int n,const int m, const HGArrayBigReal &imass, const HGArrayBigReal &length2, const HGArrayInt &ial, const HGArrayInt &ibl, const HGArrayVector &refab, const BigReal tolf, const int ntrial);

void mollify(CompAtom *qtilde,const HGArrayVector &q0,const BigReal *lambda, HGArrayVector &force,const int n, const int m, const HGArrayBigReal &imass,const HGArrayInt &ial,const HGArrayInt &ibl,const HGArrayVector &refab);


// DMK - Atom Separation (water vs. non-water)
#if NAMD_SeparateWaters != 0

// Macro to test if a hydrogen group represents a water molecule.
// NOTE: This test is the same test in Molecule.C for setting the
//   OxygenAtom flag in status.
// hgtype should be the number of atoms in a water hydrogen group
// It must now be set based on simulation parameters because we might
// be using tip4p

// DJH: This will give false positive for full Drude model,
//      e.g. O D H is not water but has hgs==3
#define IS_HYDROGEN_GROUP_WATER(hgs, mass)                 \
  ((hgs >= 3) && ((mass >= 14.0) && (mass <= 18.0)))

#endif


HomePatch::HomePatch(PatchID pd, int atomCnt) : Patch(pd)
// DMK - Atom Separation (water vs. non-water)
#if NAMD_SeparateWaters != 0
  ,tempAtom()
#endif
{
  min.x = PatchMap::Object()->min_a(patchID);
  min.y = PatchMap::Object()->min_b(patchID);
  min.z = PatchMap::Object()->min_c(patchID);
  max.x = PatchMap::Object()->max_a(patchID);
  max.y = PatchMap::Object()->max_b(patchID);
  max.z = PatchMap::Object()->max_c(patchID);
  center = 0.5*(min+max);

  aAway = PatchMap::Object()->numaway_a();
  bAway = PatchMap::Object()->numaway_b();
  cAway = PatchMap::Object()->numaway_c();

  migrationSuspended = false;
  allMigrationIn = false;
  marginViolations = 0;
  patchMapRead = 0; // We delay read of PatchMap data
		    // to make sure it is really valid
  inMigration = false;
  numMlBuf = 0;
  flags.sequence = -1;

  numAtoms = atomCnt;
  replacementForces = 0;

  SimParameters *simParams = Node::Object()->simParameters;
  doPairlistCheck_newTolerance = 
	0.5 * ( simParams->pairlistDist - simParams->cutoff );


  numFixedAtoms = 0;
  //if ( simParams->fixedAtomsOn ) {
  //  for ( int i = 0; i < numAtoms; ++i ) {
  //    numFixedAtoms += ( atom[i].atomFixed ? 1 : 0 );
  //  }
  //}

#ifdef NODEAWARE_PROXY_SPANNINGTREE
  ptnTree.resize(0);
  /*children = NULL;
  numChild = 0;*/
#else
  child =  new int[proxySpanDim];
  nChild = 0;	// number of proxy spanning tree children
#endif

#if CMK_PERSISTENT_COMM
  phsReady = 0;
  nphs = 0;
  localphs = new PersistentHandle[CkNumPes()];
  for (int i=0; i<CkNumPes(); i++) localphs[i] = 0;
#endif


  // DMK - Atom Separation (water vs. non-water)
  #if NAMD_SeparateWaters != 0

    // Create the scratch memory for separating atoms
    tempAtom.resize(numAtoms);
    numWaterAtoms = -1;

  #endif
  // Handle unusual water models here
  if (simParams->watmodel == WAT_TIP4) init_tip4();
  else if (simParams->watmodel == WAT_SWM4) init_swm4();

  isNewProxyAdded = 0;
}

HomePatch::HomePatch(PatchID pd, FullAtomList al) : Patch(pd), atom(al)
// DMK - Atom Separation (water vs. non-water)
#if NAMD_SeparateWaters != 0
  ,tempAtom()
#endif
{ 
  min.x = PatchMap::Object()->min_a(patchID);
  min.y = PatchMap::Object()->min_b(patchID);
  min.z = PatchMap::Object()->min_c(patchID);
  max.x = PatchMap::Object()->max_a(patchID);
  max.y = PatchMap::Object()->max_b(patchID);
  max.z = PatchMap::Object()->max_c(patchID);
  center = 0.5*(min+max);

  aAway = PatchMap::Object()->numaway_a();
  bAway = PatchMap::Object()->numaway_b();
  cAway = PatchMap::Object()->numaway_c();

  migrationSuspended = false;
  allMigrationIn = false;
  marginViolations = 0;
  patchMapRead = 0; // We delay read of PatchMap data
		    // to make sure it is really valid
  inMigration = false;
  numMlBuf = 0;
  flags.sequence = -1;

  numAtoms = atom.size();
  replacementForces = 0;

  SimParameters *simParams = Node::Object()->simParameters;
  doPairlistCheck_newTolerance = 
	0.5 * ( simParams->pairlistDist - simParams->cutoff );


  numFixedAtoms = 0;
  if ( simParams->fixedAtomsOn ) {
    for ( int i = 0; i < numAtoms; ++i ) {
      numFixedAtoms += ( atom[i].atomFixed ? 1 : 0 );
    }
  }

#ifdef NODEAWARE_PROXY_SPANNINGTREE
  ptnTree.resize(0);
  /*children = NULL;
  numChild = 0;*/
#else
  child =  new int[proxySpanDim];
  nChild = 0;	// number of proxy spanning tree children
#endif

#if CMK_PERSISTENT_COMM
  phsReady = 0;
  nphs = 0;
  localphs = new PersistentHandle[CkNumPes()];
  for (int i=0; i<CkNumPes(); i++) localphs[i] = 0;
#endif


  // DMK - Atom Separation (water vs. non-water)
  #if NAMD_SeparateWaters != 0

    // Create the scratch memory for separating atoms
    tempAtom.resize(numAtoms);
    numWaterAtoms = -1;

    // Separate the current list of atoms
    separateAtoms();

  #endif
    
  // Handle unusual water models here
  if (simParams->watmodel == WAT_TIP4) init_tip4();
  else if (simParams->watmodel == WAT_SWM4) init_swm4();

  isNewProxyAdded = 0;

}

void HomePatch::write_tip4_props() {
  printf("Writing r_om and r_ohc: %f | %f\n", r_om, r_ohc);
}

void HomePatch::init_tip4() {
  // initialize the distances needed for the tip4p water model

  Molecule *mol = Node::Object()->molecule;
  r_om = mol->r_om;
  r_ohc = mol->r_ohc;
}


void ::HomePatch::init_swm4() {
  // initialize the distances needed for the SWM4 water model
  SimParameters *simParams = Node::Object()->simParameters;
  Molecule *mol = Node::Object()->molecule;
  int ig;
  if (RIGID_NONE == simParams->rigidBonds) return;
  if (WAT_SWM4 != simParams->watmodel) return;
  for (ig = 0;  ig < numAtoms;  ig += atom[ig].hydrogenGroupSize ) {
    // find a water
    if (mol->rigid_bond_length(atom[ig].id) > 0) {
      // water is guaranteed by Molecule to have order:  O D LP H1 H2
      BigReal r_hh = mol->rigid_bond_length(atom[ig].id);
      BigReal r_oh = mol->rigid_bond_length(atom[ig+3].id);
      r_om = mol->rigid_bond_length(atom[ig+2].id);
      r_ohc = sqrt(r_oh * r_oh - 0.25 * r_hh * r_hh);
      //printf("r_om and r_ohc initialized to %f and %f\n", r_om, r_ohc);
      break;
    }
  }
}

void HomePatch::reinitAtoms(FullAtomList al) {
  atom = al;
  numAtoms = atom.size();

  // DMK - Atom Separation (water vs. non-water)
  #if NAMD_SeparateWaters != 0

    // Reset the numWaterAtoms value
    numWaterAtoms = -1;

    // Separate the atoms
    separateAtoms();

  #endif
}

// Bind a Sequencer to this HomePatch
void HomePatch::useSequencer(Sequencer *sequencerPtr)
{ sequencer=sequencerPtr; }

// start simulation over this Patch of atoms
void HomePatch::runSequencer(void)
{ sequencer->run(); }

void HomePatch::readPatchMap() {
  // iout << "Patch " << patchID << " has " << proxy.size() << " proxies.\n" << endi;
  PatchMap *p = PatchMap::Object();
  PatchID nnPatchID[PatchMap::MaxOneAway];

  patchMigrationCounter = numNeighbors 
    = PatchMap::Object()->oneAwayNeighbors(patchID, nnPatchID);
  DebugM( 1, "NumNeighbors for pid " <<patchID<<" is "<< numNeighbors << "\n");
  int n;
  for (n=0; n<numNeighbors; n++) {
    realInfo[n].destNodeID = p->node(realInfo[n].destPatchID = nnPatchID[n]);
     DebugM( 1, " nnPatchID=" <<nnPatchID[n]<<" nnNodeID="<< realInfo[n].destNodeID<< "\n");
    realInfo[n].mList.resize(0);
  }

  // Make mapping from the 3x3x3 cube of pointers to real migration info
  for (int i=0; i<3; i++)
    for (int j=0; j<3; j++)
      for (int k=0; k<3; k++)
      {
	int pid =  p->pid(p->index_a(patchID)+i-1, 
	    p->index_b(patchID)+j-1, p->index_c(patchID)+k-1);
	if (pid < 0) {
	   DebugM(5, "ERROR, for patchID " << patchID <<" I got neigh pid = " << pid << "\n");
	}
	if (pid == patchID && ! (
		( (i-1) && p->periodic_a() ) ||
		( (j-1) && p->periodic_b() ) ||
		( (k-1) && p->periodic_c() ) )) {
	  mInfo[i][j][k] = NULL;
	}
	else {
	  // Does not work as expected for periodic with only two patches.
	  // Also need to check which image we want, but OK for now.  -JCP
	  for (n = 0; n<numNeighbors; n++) {
	    if (pid == realInfo[n].destPatchID) {
	      mInfo[i][j][k] = &realInfo[n];
	      break;
	    }
	  }
	  if (n == numNeighbors) { // disaster! 
	    DebugM(4,"BAD News, I could not find PID " << pid << "\n");
	  }
	}
      }

  DebugM(1,"Patch("<<patchID<<") # of neighbors = " << numNeighbors << "\n");
}

HomePatch::~HomePatch()
{
#ifdef NODEAWARE_PROXY_SPANNINGTREE
    ptnTree.resize(0);
    delete [] children;
    #ifdef USE_NODEPATCHMGR
    delete [] nodeChildren;    
    #endif
#else
  delete [] child;
#endif
}


void HomePatch::boxClosed(int)
{
  if ( ! --boxesOpen )
  {
    if ( replacementForces ) {
      for ( int i = 0; i < numAtoms; ++i ) {
        if ( replacementForces[i].replace ) {
          for ( int j = 0; j < Results::maxNumForces; ++j ) { f[j][i] = 0; }
          f[Results::normal][i] = replacementForces[i].force;
        }
      }
      replacementForces = 0;
    }
    DebugM(1,patchID << ": " << CthSelf() << " awakening sequencer "
	<< sequencer->thread << "(" << patchID << ") @" << CmiTimer() << "\n");
    // only awaken suspended threads.  Then say it is suspended.
    sequencer->awaken();
    return;
  }
  else
  {
    DebugM(1,patchID << ": " << boxesOpen << " boxes left to close.\n");
  }
}

void HomePatch::registerProxy(RegisterProxyMsg *msg) {
  DebugM(4, "registerProxy("<<patchID<<") - adding node " <<msg->node<<"\n");
  proxy.add(ProxyListElem(msg->node,forceBox.checkOut()));

  isNewProxyAdded = 1;

  Random((patchID + 37) * 137).reorder(proxy.begin(),proxy.size());
  delete msg;
}

void HomePatch::unregisterProxy(UnregisterProxyMsg *msg) {
  int n = msg->node;
  ProxyListElem *pe = proxy.begin();
  for ( ; pe->node != n; ++pe );
  forceBox.checkIn(pe->forceBox);
  proxy.del(pe - proxy.begin());
  delete msg;
}

#if USE_TOPOMAP && USE_SPANNING_TREE

int HomePatch::findSubroots(int dim, int* subroots, int psize, int* pidscopy){
  int nChild = 0;
  int cones[6][proxySpanDim*proxySpanDim+proxySpanDim];
  int conesizes[6] = {0,0,0,0,0,0};
  int conecounters[6] = {0,0,0,0,0,0};
  int childcounter = 0;
  nChild = (psize>proxySpanDim)?proxySpanDim:psize;
  TopoManager tmgr;
  for(int i=0;i<psize;i++){
    int cone = tmgr.getConeNumberForRank(pidscopy[i]);
    cones[cone][conesizes[cone]++] = pidscopy[i];
  }

  while(childcounter<nChild){
    for(int i=0;i<6;i++){
      if(conecounters[i]<conesizes[i]){
        subroots[childcounter++] = cones[i][conecounters[i]++];
      }
    }
  }
  for(int i=nChild;i<proxySpanDim;i++)
    subroots[i] = -1;
  return nChild;
}
#endif // USE_TOPOMAP 

static int compDistance(const void *a, const void *b)
{
  int d1 = abs(*(int *)a - CkMyPe());
  int d2 = abs(*(int *)b - CkMyPe());
  if (d1 < d2) 
    return -1;
  else if (d1 == d2) 
    return 0;
  else 
    return 1;
}

void HomePatch::sendProxies()
{
  ProxyListIter pli(proxy);
  NodeIDList list;
  for ( pli = pli.begin(); pli != pli.end(); ++pli )
  {
    list.add(pli->node);
  }
  ProxyMgr::Object()->sendProxies(patchID, &list[0], list.size());  
}

#ifdef NODEAWARE_PROXY_SPANNINGTREE
void HomePatch::buildNodeAwareSpanningTree(void){
    //build the naive spanning tree for this home patch    
    int *proxyNodeMap = new int[CkNumNodes()]; //each element indiates the number of proxies residing on this node 
    NodeIDList proxyPeList;
    int proxyCnt = proxy.size();
    if(proxyCnt==0) {
        //this case will not happen in practice.
        //In debugging state where spanning tree is enforced, then this could happen
        //Chao Mei        
        #if defined(PROCTRACE_DEBUG) && defined(NAST_DEBUG)
        DebugFileTrace *dft = DebugFileTrace::Object();
        dft->openTrace();
        dft->writeTrace("HomePatch[%d] has 0 proxy on proc[%d] node[%d]\n", patchID, CkMyPe(), CkMyNode());
        dft->closeTrace();
        #endif
        return;
    }
    proxyPeList.resize(proxyCnt);
    ProxyListIter pli(proxy);
    int i=0;
    for ( pli = pli.begin(); pli != pli.end(); ++pli, i++ )
        proxyPeList[i] = pli->node;    
    ProxyMgr::buildSinglePatchNodeAwareSpanningTree(patchID, proxyPeList, ptnTree, proxyNodeMap);    
    delete [] proxyNodeMap;
    proxyPeList.resize(0);

    //optimize on the naive spanning tree

    //setup the children
    setupChildrenFromProxySpanningTree();
    //send down to children
    sendNodeAwareSpanningTree();
}

void HomePatch::setupChildrenFromProxySpanningTree(){
    if(ptnTree.size()==0) {
        numChild = 0;
        delete [] children;
        children = NULL;
        #ifdef USE_NODEPATCHMGR
        numNodeChild = 0;
        delete [] nodeChildren;
        nodeChildren = NULL;        
        #endif
        return;
    }
    proxyTreeNode *rootnode = &ptnTree.item(0);
    CmiAssert(rootnode->peIDs[0] == CkMyPe());
    //set up children
    //1. add external children (the first proc inside the proxy tree node)    
    //2. add internal children (with threshold that would enable spanning    
    int internalChild;
    int externalChild;
    internalChild = rootnode->numPes-1;
    numChild = internalChild;
    if(numChild > inNodeProxySpanDim) {        
        //tree construction within a node)
        CmiAbort("Enabling in-node spanning tree construction has not been implemented yet!\n");
    }else{
        //exclude the root node
        int treesize = ptnTree.size();
        externalChild = (proxySpanDim>(treesize-1))?(treesize-1):proxySpanDim;
        numChild += externalChild;        

        delete [] children;
        #ifdef USE_NODEPATCHMGR
        delete [] nodeChildren;
        #endif
        if(numChild==0){
            children = NULL;
            #ifdef USE_NODEPATCHMGR
            nodeChildren = NULL;
            numNodeChild = 0;
            #endif
            return;
        }
        children = new int[numChild];    
        for(int i=0; i<externalChild; i++) {
            children[i] = ptnTree.item(i+1).peIDs[0];
        }
        for(int i=externalChild, j=1; i<numChild; i++, j++) {
            children[i] = rootnode->peIDs[j];
        }
    }

    #ifdef USE_NODEPATCHMGR
    //only register the cores that have proxy patches. The HomePach's core
    //doesn't need to be registered.
    CProxy_NodeProxyMgr pm(CkpvAccess(BOCclass_group).nodeProxyMgr);
    NodeProxyMgr *npm = pm[CkMyNode()].ckLocalBranch();
    if(rootnode->numPes==1){
        npm->registerPatch(patchID, 0, NULL);        
    }
    else{
        npm->registerPatch(patchID, rootnode->numPes-1, &rootnode->peIDs[1]);        
    }

    //set up childrens in terms of node ids
    numNodeChild = externalChild;
    if(internalChild) numNodeChild++;
    nodeChildren = new int[numNodeChild];    
    for(int i=0; i<externalChild; i++) {
        nodeChildren[i] = ptnTree.item(i+1).nodeID;        
    }
    //the last entry always stores this node id if there are proxies
    //on other cores of the same node for this patch.
    if(internalChild)
        nodeChildren[numNodeChild-1] = rootnode->nodeID;
    #endif
    
    #if defined(PROCTRACE_DEBUG) && defined(NAST_DEBUG)
    DebugFileTrace *dft = DebugFileTrace::Object();
    dft->openTrace();
    dft->writeTrace("HomePatch[%d] has %d children: ", patchID, numChild);
    for(int i=0; i<numChild; i++)
        dft->writeTrace("%d ", children[i]);
    dft->writeTrace("\n");
    dft->closeTrace();
    #endif
    
}
#endif

#ifdef NODEAWARE_PROXY_SPANNINGTREE
//This is not an entry method, but takes an argument of message type
void HomePatch::recvNodeAwareSpanningTree(ProxyNodeAwareSpanningTreeMsg *msg){
    //set up the whole tree ptnTree
    int treesize = msg->numNodesWithProxies;    
    ptnTree.resize(treesize);    
    int *pAllPes = msg->allPes;
    for(int i=0; i<treesize; i++) {
        proxyTreeNode *oneNode = &ptnTree.item(i);
        delete oneNode->peIDs;
        oneNode->numPes = msg->numPesOfNode[i];
        oneNode->nodeID = CkNodeOf(*pAllPes);
        oneNode->peIDs = new int[oneNode->numPes];
        for(int j=0; j<oneNode->numPes; j++) {
            oneNode->peIDs[j] = *pAllPes;
            pAllPes++;
        }
    }
    //setup children
    setupChildrenFromProxySpanningTree();
    //send down to children
    sendNodeAwareSpanningTree();
}

void HomePatch::sendNodeAwareSpanningTree(){
    if(ptnTree.size()==0) return;    
    ProxyNodeAwareSpanningTreeMsg *msg = 
        ProxyNodeAwareSpanningTreeMsg::getANewMsg(patchID, CkMyPe(), ptnTree.begin(), ptnTree.size());
   
    #if defined(PROCTRACE_DEBUG) && defined(NAST_DEBUG)
    msg->printOut("HP::sendST");
    #endif

    CmiAssert(CkMyPe() == msg->allPes[0]);
    ProxyMgr::Object()->sendNodeAwareSpanningTree(msg);

}
#else
void HomePatch::recvNodeAwareSpanningTree(ProxyNodeAwareSpanningTreeMsg *msg){}
void HomePatch::sendNodeAwareSpanningTree(){}
#endif

#ifndef NODEAWARE_PROXY_SPANNINGTREE
// recv a spanning tree from load balancer
void HomePatch::recvSpanningTree(int *t, int n)
{
  int i;
  nChild=0;
  tree.resize(n);
  for (i=0; i<n; i++) {
    tree[i] = t[i];
  }

  for (i=1; i<=proxySpanDim; i++) {
    if (tree.size() <= i) break;
    child[i-1] = tree[i];
    nChild++;
  }

  // send down to children
  sendSpanningTree();
}

void HomePatch::sendSpanningTree()
{
  if (tree.size() == 0) return;
  ProxySpanningTreeMsg *msg = new ProxySpanningTreeMsg;
  msg->patch = patchID;
  msg->node = CkMyPe();
  msg->tree = tree;
  ProxyMgr::Object()->sendSpanningTree(msg);  
}
#else
void HomePatch::recvSpanningTree(int *t, int n){}
void HomePatch::sendSpanningTree(){}
#endif

#ifndef NODEAWARE_PROXY_SPANNINGTREE
void HomePatch::buildSpanningTree(void)
{
  nChild = 0;
  int psize = proxy.size();
  if (psize == 0) return;
  NodeIDList oldtree = tree;
  int oldsize = oldtree.size();
  tree.resize(psize + 1);
  tree.setall(-1);
  tree[0] = CkMyPe();
  int s=1, e=psize+1;
  ProxyListIter pli(proxy);
  int patchNodesLast =
    ( PatchMap::Object()->numNodesWithPatches() < ( 0.7 * CkNumPes() ) );
  int nNonPatch = 0;
#if 0
  int* pelists = new int[psize];
  for (int i=0; i<psize; i++) pelists[i] = -1;
  for ( pli = pli.begin(); pli != pli.end(); ++pli ) {
    int idx = rand()%psize;
    while (pelists[idx] != -1) { idx++; if (idx == psize) idx = 0; }
    pelists[idx] = pli->node;
  }
  for ( int i=0; i<psize; i++)
  {
    if ( patchNodesLast && PatchMap::Object()->numPatchesOnNode(pelists[i]) ) {
      tree[--e] = pelists[i];
    } else {
      tree[s++] = pelists[i];
      nNonPatch++;
    }
  }
  delete [] pelists;
#else
    // try to put it to the same old tree
  for ( pli = pli.begin(); pli != pli.end(); ++pli )
  {
    int oldindex = oldtree.find(pli->node);
    if (oldindex != -1 && oldindex < psize) {
      tree[oldindex] = pli->node;
    }
  }
  s=1; e=psize;
  for ( pli = pli.begin(); pli != pli.end(); ++pli )
  {
    if (tree.find(pli->node) != -1) continue;    // already assigned
    if ( patchNodesLast && PatchMap::Object()->numPatchesOnNode(pli->node) ) {
      while (tree[e] != -1) { e--; if (e==-1) e = psize; }
      tree[e] = pli->node;
    } else {
      while (tree[s] != -1) { s++; if (s==psize+1) s = 1; }
      tree[s] = pli->node;
      nNonPatch++;
    }
  }
#if 1
  if (oldsize==0 && nNonPatch) {
    // first time, sort by distance
    qsort(&tree[1], nNonPatch, sizeof(int), compDistance);
  }
#endif

  //CkPrintf("home: %d:(%d) %d %d %d %d %d\n", patchID, tree.size(),tree[0],tree[1],tree[2],tree[3],tree[4]);

#if USE_TOPOMAP && USE_SPANNING_TREE
  
  //Right now only works for spanning trees with two levels
  int *treecopy = new int [psize];
  int subroots[proxySpanDim];
  int subsizes[proxySpanDim];
  int subtrees[proxySpanDim][proxySpanDim];
  int idxes[proxySpanDim];
  int i = 0;

  for(i=0;i<proxySpanDim;i++){
    subsizes[i] = 0;
    idxes[i] = i;
  }
  
  for(i=0;i<psize;i++){
    treecopy[i] = tree[i+1];
  }
  
  TopoManager tmgr;
  tmgr.sortRanksByHops(treecopy,nNonPatch);
  tmgr.sortRanksByHops(treecopy+nNonPatch,
						psize-nNonPatch);  
  
  /* build tree and subtrees */
  nChild = findSubroots(proxySpanDim,subroots,psize,treecopy);
  delete [] treecopy;
  
  for(int i=1;i<psize+1;i++){
    int isSubroot=0;
    for(int j=0;j<nChild;j++)
      if(tree[i]==subroots[j]){
        isSubroot=1;
	break;
      }
    if(isSubroot) continue;
    
    int bAdded = 0;
    tmgr.sortIndexByHops(tree[i], subroots,
						  idxes, proxySpanDim);
    for(int j=0;j<proxySpanDim;j++){
      if(subsizes[idxes[j]]<proxySpanDim){
        subtrees[idxes[j]][(subsizes[idxes[j]])++] = tree[i];
	bAdded = 1; 
        break;
      }
    }
    if( psize > proxySpanDim && ! bAdded ) {
      NAMD_bug("HomePatch BGL Spanning Tree error: Couldn't find subtree for leaf\n");
    }
  }

#else /* USE_TOPOMAP && USE_SPANNING_TREE */
  
  for (int i=1; i<=proxySpanDim; i++) {
    if (tree.size() <= i) break;
    child[i-1] = tree[i];
    nChild++;
  }
#endif
#endif
  
#if 0
  // for debugging
  CkPrintf("[%d] Spanning tree for %d with %d children %d nNonPatch %d\n", CkMyPe(), patchID, psize, nNonPatch);
  for (int i=0; i<psize+1; i++) {
    CkPrintf("%d ", tree[i]);
  }
  CkPrintf("\n");
#endif
  // send to children nodes
  sendSpanningTree();
}
#endif


void HomePatch::receiveResults(ProxyResultVarsizeMsg *msg){
    DebugM(4, "patchID("<<patchID<<") receiveRes() nodeID("<<msg->node<<")\n");
    int n = msg->node;
    ProxyListElem *pe = proxy.begin();
    for ( ; pe->node != n; ++pe );
    Results *r = pe->forceBox->open();

    char *iszeroPtr = msg->isZero;
    Force *msgFPtr = msg->forceArr;

    for ( int k = 0; k < Results::maxNumForces; ++k )
    {
      Force *rfPtr = r->f[k];
      for(int i=0; i<msg->flLen[k]; i++, rfPtr++, iszeroPtr++) {
          if((*iszeroPtr)!=1) {
              *rfPtr += *msgFPtr;
              msgFPtr++;
          }
      }      
    }
    pe->forceBox->close(&r);
    delete msg;
}

void HomePatch::receiveResults(ProxyResultMsg *msg)
{
  DebugM(4, "patchID("<<patchID<<") receiveRes() nodeID("<<msg->node<<")\n");
  int n = msg->node;
  ProxyListElem *pe = proxy.begin();
  for ( ; pe->node != n; ++pe );
  Results *r = pe->forceBox->open();
  for ( int k = 0; k < Results::maxNumForces; ++k )
  {
    Force *f = r->f[k];
    register ForceList::iterator f_i, f_e;
    f_i = msg->forceList[k].begin();
    f_e = msg->forceList[k].end();
    for ( ; f_i != f_e; ++f_i, ++f ) *f += *f_i;
  }
  pe->forceBox->close(&r);
  delete msg;
}

void HomePatch::receiveResults(ProxyCombinedResultMsg *msg)
{
  DebugM(4, "patchID("<<patchID<<") receiveRes() #nodes("<<msg->nodes.size()<<")\n");
//CkPrintf("[%d] Homepatch: %d receiveResults from %d nodes\n", CkMyPe(), patchID, n);
  for (int i=0; i<msg->nodes.size(); i++) {
    int node = msg->nodes[i];
    ProxyListElem *pe = proxy.begin();
    for ( ; pe->node != node; ++pe );
    Results *r = pe->forceBox->open();
    if (i==0) {
      for ( int k = 0; k < Results::maxNumForces; ++k )
      {
        Force *f = r->f[k];
        register ForceList::iterator f_i, f_e;
        f_i = msg->forceList[k].begin();
        f_e = msg->forceList[k].end();
        //for ( ; f_i != f_e; ++f_i, ++f ) *f += *f_i;

	int nf = f_e - f_i;
#ifdef ARCH_POWERPC
#pragma disjoint (*f_i, *f)
#pragma unroll(4)
#endif
	for (int count = 0; count < nf; count++) {
	  f[count].x += f_i[count].x;      
	  f[count].y += f_i[count].y;      
	  f[count].z += f_i[count].z;
	}
      }
    }
    pe->forceBox->close(&r);
  }

    delete msg;
}

void HomePatch::positionsReady(int doMigration)
{    
  flags.sequence++;

  if (!patchMapRead) { readPatchMap(); }
      
  if (numNeighbors) {
    if (doMigration) {
      doAtomMigration();
    } else {
      doMarginCheck();
    }
  }
  doMigration = (doMigration && numNeighbors) || ! patchMapRead;

  // Workaround for oversize groups
  doGroupSizeCheck();

  // Copy information needed by computes and proxys to Patch::p.
  p.resize(numAtoms);
  CompAtom *p_i = p.begin();
  pExt.resize(numAtoms);
  CompAtomExt *pExt_i = pExt.begin();
  FullAtom *a_i = atom.begin();
  int i; int n = numAtoms;
  for ( i=0; i<n; ++i ) { 
    p_i[i] = a_i[i]; 
    pExt_i[i] = a_i[i];
  }

  // Measure atom movement to test pairlist validity
  doPairlistCheck();

  if (flags.doMolly) mollyAverage();

  // Must Add Proxy Changes when migration completed!
  ProxyListIter pli(proxy);
  int *pids = NULL;
  int npid;
  if (proxySendSpanning == 0) {
    npid = proxy.size();
    pids = new int[npid];
    int *pidi = pids;
    int *pide = pids + proxy.size();
    int patchNodesLast =
      ( PatchMap::Object()->numNodesWithPatches() < ( 0.7 * CkNumPes() ) );
    for ( pli = pli.begin(); pli != pli.end(); ++pli )
    {
      if ( patchNodesLast && PatchMap::Object()->numPatchesOnNode(pli->node) ) {
        *(--pide) = pli->node;
      } else {
        *(pidi++) = pli->node;
      }
    }
  }
  else {
#ifdef NODEAWARE_PROXY_SPANNINGTREE
    #ifdef USE_NODEPATCHMGR
    npid = numNodeChild;
    if(numNodeChild>0) {
        pids = new int[npid];
        memcpy(pids, nodeChildren, sizeof(int)*npid);
    }
    #else
    npid = numChild;
    if(numChild>0) {
        pids = new int[numChild];
        memcpy(pids, children, sizeof(int)*numChild);
    }
    #endif
#else
    npid = nChild;
    pids = new int[proxySpanDim];
    for (int i=0; i<nChild; i++) pids[i] = child[i];
#endif
  }
  if (npid) {
#if CMK_PERSISTENT_COMM
    if (phsReady == 0)
      {
//CmiPrintf("Build on %d phs0:%d\n", CkMyPe(), localphs[0]);
     for (int i=0; i<npid; i++) {
       localphs[i] = CmiCreatePersistent(pids[i], 30000);
     }
     nphs = npid;
     phsReady = 1;
    }
#endif
    int seq = flags.sequence;
    int priority = PROXY_DATA_PRIORITY + PATCH_PRIORITY(patchID);
    //begin to prepare proxy msg and send it
    int pdMsgPLLen = p.size();
    int pdMsgAvgPLLen = 0;
    if(flags.doMolly) {
        pdMsgAvgPLLen = p_avg.size();
    }
    int pdMsgPLExtLen = 0;
    if(doMigration || isNewProxyAdded) {
        pdMsgPLExtLen = pExt.size();
    }
    ProxyDataMsg *nmsg = new (pdMsgPLLen, pdMsgAvgPLLen, pdMsgPLExtLen, PRIORITY_SIZE) ProxyDataMsg;
    SET_PRIORITY(nmsg,seq,priority);
    nmsg->patch = patchID;
    nmsg->flags = flags;
    nmsg->plLen = pdMsgPLLen;                
    //copying data to the newly created msg
    memcpy(nmsg->positionList, p.begin(), sizeof(CompAtom)*pdMsgPLLen);
    nmsg->avgPlLen = pdMsgAvgPLLen;        
    if(flags.doMolly) {
        memcpy(nmsg->avgPositionList, p_avg.begin(), sizeof(CompAtom)*pdMsgAvgPLLen);
    }
    nmsg->plExtLen = pdMsgPLExtLen;
    if(doMigration || isNewProxyAdded){     
        memcpy(nmsg->positionExtList, pExt.begin(), sizeof(CompAtomExt)*pdMsgPLExtLen);
    }
    
#if NAMD_SeparateWaters != 0
    //DMK - Atom Separation (water vs. non-water)
    nmsg->numWaterAtoms = numWaterAtoms;
#endif

#if defined(NODEAWARE_PROXY_SPANNINGTREE) && defined(USE_NODEPATCHMGR) && (CMK_SMP) && defined(NAMDSRC_IMMQD_HACK)
    nmsg->isFromImmMsgCall = 0;
#endif
    
    #if defined(PROCTRACE_DEBUG) && defined(NAST_DEBUG)
    DebugFileTrace *dft = DebugFileTrace::Object();
    dft->openTrace();
    dft->writeTrace("HP::posReady: for HomePatch[%d], sending proxy msg to: ", patchID);
    for(int i=0; i<npid; i++) {
        dft->writeTrace("%d ", pids[i]);
    }
    dft->writeTrace("\n");
    dft->closeTrace();
    #endif

    if(doMigration) {
        ProxyMgr::Object()->sendProxyAll(nmsg,npid,pids);
    }else{
        ProxyMgr::Object()->sendProxyData(nmsg,npid,pids);
    }
#if CMK_PERSISTENT_COMM
    CmiUsePersistentHandle(NULL, 0);
#endif
    isNewProxyAdded = 0;
  }
  delete [] pids;
  DebugM(4, "patchID("<<patchID<<") doing positions Ready\n");

#ifdef REMOVE_PROXYDATAMSG_EXTRACOPY
  positionPtrBegin = p.begin();
  positionPtrEnd = p.end();
#endif

  if(flags.doMolly) {
      avgPositionPtrBegin = p_avg.begin();
      avgPositionPtrEnd = p_avg.end();
  }
  Patch::positionsReady(doMigration);

  patchMapRead = 1;

  // gzheng
  if (useSync) Sync::Object()->PatchReady();
}

void HomePatch::replaceForces(ExtForce *f)
{
  replacementForces = f;
}

void HomePatch::saveForce(const int ftag)
{
  f_saved[ftag].resize(numAtoms);
  for ( int i = 0; i < numAtoms; ++i )
  {
    f_saved[ftag][i] = f[ftag][i];
  }
}

#undef DEBUG_REDISTRIB_FORCE 
#undef DEBUG_REDISTRIB_FORCE_VERBOSE
/* Redistribute forces from the massless lonepair charge particle onto
 * the other atoms of the water.
 *
 * This is done using the same algorithm as charmm uses for TIP4P lonepairs.
 *
 * Pass by reference the forces (O H1 H2 LP) to be modified,
 * pass by constant reference the corresponding positions,
 * and a pointer to virial.
 */
void HomePatch::redistrib_lp_force(
    Vector& f_ox, Vector& f_h1, Vector& f_h2, Vector& f_lp,
    const Vector& p_ox, const Vector& p_h1, const Vector& p_h2,
    const Vector& p_lp, Tensor *virial) {

  Tensor wc;  // accumulate virial contribution from force redistribution

#ifdef DEBUG_REDISTRIB_FORCE 
  // Debug information to check against results at end

  // total force and torque relative to origin
  Vector totforce(0.0, 0.0, 0.0);
  Vector tottorque(0.0, 0.0, 0.0);

  totforce = f_ox + f_h1 + f_h2 + f_lp;

  tottorque += cross(f_ox, p_ox);
  tottorque += cross(f_h1, p_h1);
  tottorque += cross(f_h2, p_h2);
  //printf("Torque without LP is %f/%f/%f\n",
  //    tottorque.x, tottorque.y, tottorque.z);
  tottorque += cross(f_lp, p_lp);
  //printf("Torque with LP is %f/%f/%f\n",
  //    tottorque.x, tottorque.y, tottorque.z);
#endif

  // Calculate the radial component of the force and add it to the oxygen
  Vector r_ox_lp = p_lp - p_ox;
  BigReal rad_factor = (f_lp * r_ox_lp) * r_ox_lp.rlength() * r_ox_lp.rlength();
  Vector f_rad = r_ox_lp * rad_factor;

  Tensor vir = outer(f_rad, p_ox);
  wc += vir;

  f_ox = f_ox + f_rad;

  // Calculate the angular component
  Vector r_hcom_ox = p_ox - ( (p_h1 + p_h2) * 0.5 );
  Vector r_h2_h1_2 = (p_h1 - p_h2) * 0.5; // half of r_h2_h1

  // deviation from collinearity of charge site
  Vector r_oop = cross(r_ox_lp, r_hcom_ox);
  //
  // vector out of o-h-h plane
  Vector r_perp = cross(r_hcom_ox, r_h2_h1_2);

  // Here we assume that Ox/Lp/Hcom are linear
  // If you want to correct for deviations, this is the place

//  printf("Deviation from linearity for ox %i: %f/%f/%f\n", oxind, r_oop.x, r_oop.y, r_oop.z);

  Vector f_ang = f_lp - f_rad; // leave the angular component

  // now split this component onto the other atoms
  BigReal oxcomp = (r_hcom_ox.length() - r_ox_lp.length()) *
    r_hcom_ox.rlength();
  BigReal hydcomp = 0.5 * r_ox_lp.length() * r_hcom_ox.rlength();

  f_ox = f_ox + (f_ang * oxcomp);
  f_h1 = f_h1 + (f_ang * hydcomp);
  f_h2 = f_h2 + (f_ang * hydcomp);

  // Add virial contributions
  vir = outer(f_ang * oxcomp, p_ox);
  wc += vir;
  vir = outer(f_ang * hydcomp, p_h1);
  wc += vir;
  vir = outer(f_ang * hydcomp, p_h2);
  wc += vir;
  vir = outer(-1.0 * f_lp, p_lp);
  wc += vir;

  if ( virial ) *virial += wc;

  //Vector zerovec(0.0, 0.0, 0.0);
  f_lp = Vector(0.0, 0.0, 0.0);

#ifdef DEBUG_REDISTRIB_FORCE 
  // Check that the total force and torque come out right
  Vector newforce(0.0, 0.0, 0.0);
  Vector newtorque(0.0, 0.0, 0.0);
  BigReal error = 0.0;

  newforce = f_ox + f_h1 + f_h2;

  newtorque += cross(f_ox, p_ox);
  newtorque += cross(f_h1, p_h1);
  newtorque += cross(f_h2, p_h2);

  error = fabs(newforce.length() - totforce.length());
  if (error > 0.0001) {
     printf("Error:  Force redistribution for water "
         "exceeded force tolerance (%f vs. %f)\n",
         newforce.length(), totforce.length());
  }
#ifdef DEBUG_REDISTRIB_FORCE_VERBOSE
  printf("Error in force length:  %f\n", error);
#endif

  error = fabs(newtorque.length() - tottorque.length());
  if (error > 0.0001) {
     printf("Error:  Force redistribution for water "
         "exceeded torque tolerance (%f vs. %f)\n",
         newtorque.length(), tottorque.length());
  }
#ifdef DEBUG_REDISTRIB_FORCE_VERBOSE
  printf("Error in torque:  %f\n", error);
#endif
#endif /* DEBUG */
}

void HomePatch::swm4_omrepos(Vector *ref, Vector *pos, Vector *vel,
    BigReal invdt) {
  // Reposition lonepair (Om) particle of Drude SWM4 water.
  // Same comments apply as to tip4_omrepos(), but the ordering of atoms
  // is different: O, D, LP, H1, H2.
  pos[2] = pos[0] + (0.5 * (pos[3] + pos[4]) - pos[0]) * (r_om / r_ohc);
  // Now, adjust velocity of particle to get it to appropriate place
  if (invdt != 0) {
    vel[2] = (pos[2] - ref[2]) * invdt;
  }
  // No virial correction needed since lonepair is massless
}

void HomePatch::tip4_omrepos(Vector* ref, Vector* pos, Vector* vel, BigReal invdt) {
  /* Reposition the om particle of a tip4p water
   * A little geometry shows that the appropriate position is given by
   * R_O + (1 / 2 r_ohc) * ( 0.5 (R_H1 + R_H2) - R_O ) 
   * Here r_om is the distance from the oxygen to Om site, and r_ohc
   * is the altitude from the oxygen to the hydrogen center of mass
   * Those quantities are precalculated upon initialization of HomePatch
   *
   * Ordering of TIP4P atoms: O, H1, H2, LP.
   */

  //printf("rom/rohc are %f %f and invdt is %f\n", r_om, r_ohc, invdt);
  //printf("Other positions are: \n  0: %f %f %f\n  1: %f %f %f\n  2: %f %f %f\n", pos[0].x, pos[0].y, pos[0].z, pos[1].x, pos[1].y, pos[1].z, pos[2].x, pos[2].y, pos[2].z);
  pos[3] = pos[0] + (0.5 * (pos[1] + pos[2]) - pos[0]) * (r_om / r_ohc); 
  //printf("New position for lp is %f %f %f\n", pos[3].x, pos[3].y, pos[3].z);

  // Now, adjust the velocity of the particle to get it to the appropriate place
  if (invdt != 0) {
    vel[3] = (pos[3] - ref[3]) * invdt;
  }

  // No virial correction needed, since this is a massless particle
  return;
}

void HomePatch::redistrib_swm4_forces(const int ftag, Tensor *virial) {
  // Loop over the patch's atoms and apply the appropriate corrections
  // to get all forces off of lone pairs
  ForceList *f_mod = f;
  for (int i = 0;  i < numAtoms;  i++) {
    if (atom[i].mass < 0.01) {
      // found lonepair
      redistrib_lp_force(f_mod[ftag][i-2], f_mod[ftag][i+1],
          f_mod[ftag][i+2], f_mod[ftag][i],
          atom[i-2].position, atom[i+1].position,
          atom[i+2].position, atom[i].position, virial);
    }
  }
}

void HomePatch::redistrib_tip4p_forces(const int ftag, Tensor* virial) {
  // Loop over the patch's atoms and apply the appropriate corrections
  // to get all forces off of lone pairs
  // Atom ordering:  O H1 H2 LP

  ForceList *f_mod =f;
  for (int i=0; i<numAtoms; i++) {
    if (atom[i].mass < 0.01) {
      // found lonepair
      redistrib_lp_force(f_mod[ftag][i-3], f_mod[ftag][i-2],
          f_mod[ftag][i-1], f_mod[ftag][i],
          atom[i-3].position, atom[i-2].position,
          atom[i-1].position, atom[i].position, virial);
    }
  }
}


void HomePatch::addForceToMomentum(const BigReal timestep, const int ftag,
							const int useSaved)
{
  SimParameters *simParams = Node::Object()->simParameters;
  const BigReal dt = timestep / TIMEFACTOR;
  ForceList *f_use = (useSaved ? f_saved : f);

  if ( simParams->fixedAtomsOn ) {
    for ( int i = 0; i < numAtoms; ++i ) {
      if ( atom[i].atomFixed ) {
        atom[i].velocity = 0;
      } else {
        BigReal recip_val = ( atom[i].mass > 0. ? dt * namd_reciprocal( atom[i].mass ) : 0.); 
        atom[i].velocity += f_use[ftag][i] * recip_val;
      }
    }
  } else {
    FullAtom *atom_arr  = atom.begin();
    const Force    *force_arr = f_use[ftag].const_begin();
#ifdef ARCH_POWERPC
#pragma disjoint (*force_arr, *atom_arr)
#endif
    for ( int i = 0; i < numAtoms; ++i ) {
      if (atom[i].mass == 0.) continue;
      BigReal recip_val = ( atom[i].mass > 0. ? dt * namd_reciprocal( atom[i].mass ) : 0.); 
      //printf("Taking reciprocal of mass %f\n", atom[i].mass);
      atom_arr[i].velocity.x += force_arr[i].x * recip_val;
      atom_arr[i].velocity.y += force_arr[i].y * recip_val;
      atom_arr[i].velocity.z += force_arr[i].z * recip_val;
    }
  }
}

void HomePatch::addVelocityToPosition(const BigReal timestep)
{
  SimParameters *simParams = Node::Object()->simParameters;
  const BigReal dt = timestep / TIMEFACTOR;
  if ( simParams->fixedAtomsOn ) {
    for ( int i = 0; i < numAtoms; ++i ) {
      if ( ! atom[i].atomFixed ) atom[i].position += atom[i].velocity * dt;
    }
  } else {
    FullAtom *atom_arr  = atom.begin();
    for ( int i = 0; i < numAtoms; ++i ) {
      atom_arr[i].position.x  +=  atom_arr[i].velocity.x * dt;
      atom_arr[i].position.y  +=  atom_arr[i].velocity.y * dt;
      atom_arr[i].position.z  +=  atom_arr[i].velocity.z * dt;
    }
  }
}

//  RATTLE algorithm from Allen & Tildesley
int HomePatch::rattle1(const BigReal timestep, Tensor *virial, 
    SubmitReduction *ppreduction)
{
  Molecule *mol = Node::Object()->molecule;
  SimParameters *simParams = Node::Object()->simParameters;
  const int fixedAtomsOn = simParams->fixedAtomsOn;
  const int useSettle = simParams->useSettle;
  const BigReal dt = timestep / TIMEFACTOR;
  const BigReal invdt = (dt == 0.) ? 0. : 1.0 / dt; // precalc 1/dt
  BigReal tol2 = 2.0 * simParams->rigidTol;
  int maxiter = simParams->rigidIter;
  int dieOnError = simParams->rigidDie;
  int i, iter;
  BigReal dsq[10], tmp;
  int ial[10], ibl[10];
  Vector ref[10];  // reference position
  Vector refab[10];  // reference vector
  Vector pos[10];  // new position
  Vector vel[10];  // new velocity
  Vector netdp[10];  // total momentum change from constraint
  BigReal rmass[10];  // 1 / mass
  int fixed[10];  // is atom fixed?
  Tensor wc;  // constraint virial
  BigReal idz, zmin;
  int nslabs;

  // Initialize the settle algorithm with water parameters
  // settle1() assumes all waters are identical,
  // and will generate bad results if they are not.
  // XXX this will move to Molecule::build_atom_status when that 
  // version is debugged
  if (!settle1isinitted()) {
    for ( int ig = 0; ig < numAtoms; ig += atom[ig].hydrogenGroupSize ) {
      // find a water
      if (mol->rigid_bond_length(atom[ig].id) > 0) {
        int oatm;
        if (simParams->watmodel == WAT_SWM4) {
          oatm = ig+3;  // skip over Drude and Lonepair
          //printf("ig=%d  mass_ig=%g  oatm=%d  mass_oatm=%g\n",
          //    ig, atom[ig].mass, oatm, atom[oatm].mass);
        }
        else {
          oatm = ig+1;
          // Avoid using the Om site to set this by mistake
          if (atom[ig].mass < 0.5 || atom[ig+1].mass < 0.5) {
            oatm += 1;
          }
        }

        // initialize settle water parameters
        settle1init(atom[ig].mass, atom[oatm].mass, 
                    mol->rigid_bond_length(atom[ig].id), 
                    mol->rigid_bond_length(atom[oatm].id));
        break; // done with init
      }
    }
  }

  if (ppreduction) {
    nslabs = simParams->pressureProfileSlabs;
    idz = nslabs/lattice.c().z;
    zmin = lattice.origin().z - 0.5*lattice.c().z;
  }

  // Size of a hydrogen group for water
  int wathgsize = 3;
  int watmodel = simParams->watmodel;
  if (watmodel == WAT_TIP4) wathgsize = 4;
  else if (watmodel == WAT_SWM4) wathgsize = 5;
  
  for ( int ig = 0; ig < numAtoms; ig += atom[ig].hydrogenGroupSize ) {
    int hgs = atom[ig].hydrogenGroupSize;
    if ( hgs == 1 ) continue;  // only one atom in group
    // cache data in local arrays and integrate positions normally
    for ( i = 0; i < hgs; ++i ) {
      ref[i] = atom[ig+i].position;
      pos[i] = atom[ig+i].position;
      vel[i] = atom[ig+i].velocity;
      rmass[i] = (atom[ig+1].mass > 0. ? 1. / atom[ig+i].mass : 0.);
      //printf("rmass of %i is %f\n", ig+i, rmass[i]);
      fixed[i] = ( fixedAtomsOn && atom[ig+i].atomFixed );
      //printf("fixed status of %i is %i\n", i, fixed[i]);
      // undo addVelocityToPosition to get proper reference coordinates
      if ( fixed[i] ) rmass[i] = 0.; else pos[i] += vel[i] * dt;
    }
    int icnt = 0;
    if ( ( tmp = mol->rigid_bond_length(atom[ig].id) ) > 0 ) {  // for water
      if (hgs != wathgsize) {
        NAMD_bug("Hydrogen group error caught in rattle1().");
      }
      // Use SETTLE for water unless some of the water atoms are fixed,
      // for speed we test groupFixed rather than the individual atoms
      if (useSettle && !atom[ig].groupFixed) {
        if (simParams->watmodel == WAT_SWM4) {
          // SWM4 ordering:  O D LP H1 H2
          // do swap(O,LP) and call settle with subarray O H1 H2
          // swap back after we return
          Vector lp_ref = ref[2];
          Vector lp_pos = pos[2];
          Vector lp_vel = vel[2];
          ref[2] = ref[0];
          pos[2] = pos[0];
          vel[2] = vel[0];
          settle1(ref+2, pos+2, vel+2, invdt);
          ref[0] = ref[2];
          pos[0] = pos[2];
          vel[0] = vel[2];
          ref[2] = lp_ref;
          pos[2] = lp_pos;
          vel[2] = lp_vel;
          // determine for LP updated pos and vel
          swm4_omrepos(ref, pos, vel, invdt);
        }
        else {
          settle1(ref, pos, vel, invdt);
          if (simParams->watmodel == WAT_TIP4) {
            tip4_omrepos(ref, pos, vel, invdt);
          }
        }

        // which slab the hydrogen group will belong to
        // for pprofile calculations.
        int ppoffset, partition;
        if ( invdt == 0 ) for ( i = 0; i < wathgsize; ++i ) {
          atom[ig+i].position = pos[i];
        } else if ( virial == 0 ) for ( i = 0; i < wathgsize; ++i ) {
          atom[ig+i].velocity = vel[i];
        } else for ( i = 0; i < wathgsize; ++i ) {
          Force df = (vel[i] - atom[ig+i].velocity) * ( atom[ig+i].mass * invdt );
          Tensor vir = outer(df, ref[i]);
          wc += vir;
          f[Results::normal][ig+i] += df;
          atom[ig+i].velocity = vel[i];
          if (ppreduction) {
            // put all the atoms from a water in the same slab.  Atom 0
            // should be the parent atom.
            if (!i) {
              BigReal z = pos[i].z;
              partition = atom[ig].partition;
              int slab = (int)floor((z-zmin)*idz);
              if (slab < 0) slab += nslabs;
              else if (slab >= nslabs) slab -= nslabs;
              ppoffset = 3*(slab + nslabs*partition);
            }
            ppreduction->item(ppoffset  ) += vir.xx;
            ppreduction->item(ppoffset+1) += vir.yy;
            ppreduction->item(ppoffset+2) += vir.zz;
          }
        }
        continue;
      }
      if ( !(fixed[1] && fixed[2]) ) {
	dsq[icnt] = tmp * tmp;  ial[icnt] = 1;  ibl[icnt] = 2;  ++icnt;
      }
    }
    for ( i = 1; i < hgs; ++i ) {  // normal bonds to mother atom
      if ( ( tmp = mol->rigid_bond_length(atom[ig+i].id) ) > 0 ) {
	if ( !(fixed[0] && fixed[i]) ) {
	  dsq[icnt] = tmp * tmp;  ial[icnt] = 0;  ibl[icnt] = i;  ++icnt;
	}
      }
    }
    if ( icnt == 0 ) continue;  // no constraints
    for ( i = 0; i < icnt; ++i ) {
      refab[i] = ref[ial[i]] - ref[ibl[i]];
    }
    for ( i = 0; i < hgs; ++i ) {
      netdp[i] = 0.;
    }
    int done;
    int consFailure;
    for ( iter = 0; iter < maxiter; ++iter ) {
//if (iter > 0) CkPrintf("iteration %d\n", iter);
      done = 1;
      consFailure = 0;
      for ( i = 0; i < icnt; ++i ) {
	int a = ial[i];  int b = ibl[i];
	Vector pab = pos[a] - pos[b];
	BigReal pabsq = pab.x*pab.x + pab.y*pab.y + pab.z*pab.z;
	BigReal rabsq = dsq[i];
	BigReal diffsq = rabsq - pabsq;
	if ( fabs(diffsq) > (rabsq * tol2) ) {
	  Vector &rab = refab[i];
	  BigReal rpab = rab.x*pab.x + rab.y*pab.y + rab.z*pab.z;
	  if ( rpab < ( rabsq * 1.0e-6 ) ) {
	    done = 0;
	    consFailure = 1;
	    continue;
	  }
	  BigReal rma = rmass[a];
	  BigReal rmb = rmass[b];
	  BigReal gab = diffsq / ( 2.0 * ( rma + rmb ) * rpab );
	  Vector dp = rab * gab;
	  pos[a] += rma * dp;
	  pos[b] -= rmb * dp;
	  if ( invdt != 0. ) {
	    dp *= invdt;
	    netdp[a] += dp;
	    netdp[b] -= dp;
	  }
	  done = 0;
	}
      }
      if ( done ) break;
    }

    if ( consFailure ) {
      if ( dieOnError ) {
	iout << iERROR << "Constraint failure in RATTLE algorithm for atom "
			<< (atom[ig].id + 1) << "!\n" << endi;
	return -1;  // triggers early exit
      } else {
	iout << iWARN << "Constraint failure in RATTLE algorithm for atom "
			<< (atom[ig].id + 1) << "!\n" << endi;
      }
    } else if ( ! done ) {
      if ( dieOnError ) {
	iout << iERROR << "Exceeded RATTLE iteration limit for atom "
			<< (atom[ig].id + 1) << "!\n" << endi;
	return -1;  // triggers early exit
      } else {
	iout << iWARN << "Exceeded RATTLE iteration limit for atom "
			<< (atom[ig].id + 1) << "!\n" << endi;
      }
    }

    // store data back to patch
    int ppoffset, partition;
    if ( invdt == 0 ) for ( i = 0; i < hgs; ++i ) {
      atom[ig+i].position = pos[i];
    } else if ( virial == 0 ) for ( i = 0; i < hgs; ++i ) {
      atom[ig+i].velocity = vel[i] + rmass[i] * netdp[i];
    } else for ( i = 0; i < hgs; ++i ) {
      Force df = netdp[i] * invdt;
      Tensor vir = outer(df, ref[i]);
      wc += vir;
      f[Results::normal][ig+i] += df;
      atom[ig+i].velocity = vel[i] + rmass[i] * netdp[i];
      if (ppreduction) {
        if (!i) {
          BigReal z = pos[i].z;
          int partition = atom[ig].partition;
          int slab = (int)floor((z-zmin)*idz);
          if (slab < 0) slab += nslabs;
          else if (slab >= nslabs) slab -= nslabs;
          ppoffset = 3*(slab + nslabs*partition);
        }
        ppreduction->item(ppoffset  ) += vir.xx;
        ppreduction->item(ppoffset+1) += vir.yy;
        ppreduction->item(ppoffset+2) += vir.zz;
      }
    }
  }
  if ( dt && virial ) *virial += wc;

  return 0;
}

//  RATTLE algorithm from Allen & Tildesley
void HomePatch::rattle2(const BigReal timestep, Tensor *virial)
{
  Molecule *mol = Node::Object()->molecule;
  SimParameters *simParams = Node::Object()->simParameters;
  const int fixedAtomsOn = simParams->fixedAtomsOn;
  const int useSettle = simParams->useSettle;
  const BigReal dt = timestep / TIMEFACTOR;
  Tensor wc;  // constraint virial
  BigReal tol = simParams->rigidTol;
  int maxiter = simParams->rigidIter;
  int dieOnError = simParams->rigidDie;
  int i, iter;
  BigReal dsqi[10], tmp;
  int ial[10], ibl[10];
  Vector ref[10];  // reference position
  Vector refab[10];  // reference vector
  Vector vel[10];  // new velocity
  BigReal rmass[10];  // 1 / mass
  BigReal redmass[10];  // reduced mass
  int fixed[10];  // is atom fixed?
  
  // Size of a hydrogen group for water
  int wathgsize = 3;
  if (simParams->watmodel == WAT_TIP4) wathgsize = 4;
  else if (simParams->watmodel == WAT_SWM4) wathgsize = 5;

  //  CkPrintf("In rattle2!\n");
  for ( int ig = 0; ig < numAtoms; ig += atom[ig].hydrogenGroupSize ) {
    //    CkPrintf("ig=%d\n",ig);
    int hgs = atom[ig].hydrogenGroupSize;
    if ( hgs == 1 ) continue;  // only one atom in group
    // cache data in local arrays and integrate positions normally
    for ( i = 0; i < hgs; ++i ) {
      ref[i] = atom[ig+i].position;
      vel[i] = atom[ig+i].velocity;
      rmass[i] = atom[ig+i].mass > 0. ? 1. / atom[ig+i].mass : 0.;
      fixed[i] = ( fixedAtomsOn && atom[ig+i].atomFixed );
      if ( fixed[i] ) rmass[i] = 0.;
    }
    int icnt = 0;
    if ( ( tmp = mol->rigid_bond_length(atom[ig].id) ) > 0 ) {  // for water
      if ( wathgsize != 4 ) {
        NAMD_bug("Hydrogen group error caught in rattle2().");
      }
      // Use SETTLE for water unless some of the water atoms are fixed
      if (useSettle && !fixed[0] && !fixed[1] && !fixed[2]) {
        settle2(atom[ig].mass, atom[ig+1].mass, ref, vel, dt, virial);
        for (i=0; i<3; i++) {
          atom[ig+i].velocity = vel[i];
        }
        continue;
      }
      if ( !(fixed[1] && fixed[2]) ) {
	redmass[icnt] = 1. / (rmass[1] + rmass[2]);
	dsqi[icnt] = 1. / (tmp * tmp);  ial[icnt] = 1;  ibl[icnt] = 2;  ++icnt;
      }
    }
    //    CkPrintf("Loop 2\n");
    for ( i = 1; i < hgs; ++i ) {  // normal bonds to mother atom
      if ( ( tmp = mol->rigid_bond_length(atom[ig+i].id) ) > 0 ) {
        if ( !(fixed[0] && fixed[i]) ) {
	  redmass[icnt] = 1. / (rmass[0] + rmass[i]);
	  dsqi[icnt] = 1. / (tmp * tmp);  ial[icnt] = 0;
	  ibl[icnt] = i;  ++icnt;
	}
      }
    }
    if ( icnt == 0 ) continue;  // no constraints
    //    CkPrintf("Loop 3\n");
    for ( i = 0; i < icnt; ++i ) {
      refab[i] = ref[ial[i]] - ref[ibl[i]];
    }
    //    CkPrintf("Loop 4\n");
    int done;
    for ( iter = 0; iter < maxiter; ++iter ) {
      done = 1;
      for ( i = 0; i < icnt; ++i ) {
	int a = ial[i];  int b = ibl[i];
	Vector vab = vel[a] - vel[b];
	Vector &rab = refab[i];
	BigReal rabsqi = dsqi[i];
	BigReal rvab = rab.x*vab.x + rab.y*vab.y + rab.z*vab.z;
	if ( (fabs(rvab) * dt * rabsqi) > tol ) {
	  Vector dp = rab * (-rvab * redmass[i] * rabsqi);
	  wc += outer(dp,rab);
	  vel[a] += rmass[a] * dp;
	  vel[b] -= rmass[b] * dp;
	  done = 0;
	}
      }
      if ( done ) break;
      //if (done) { if (iter > 0) CkPrintf("iter=%d\n", iter); break; }
    }
    if ( ! done ) {
      if ( dieOnError ) {
	NAMD_die("Exceeded maximum number of iterations in rattle2().");
      } else {
	iout << iWARN <<
	  "Exceeded maximum number of iterations in rattle2().\n" << endi;
      }
    }
    // store data back to patch
    for ( i = 0; i < hgs; ++i ) {
      atom[ig+i].velocity = vel[i];
    }
  }
  //  CkPrintf("Leaving rattle2!\n");
  // check that there isn't a constant needed here!
  *virial += wc / ( 0.5 * dt );

}


//  MOLLY algorithm part 1
void HomePatch::mollyAverage()
{
  Molecule *mol = Node::Object()->molecule;
  SimParameters *simParams = Node::Object()->simParameters;
  BigReal tol = simParams->mollyTol;
  int maxiter = simParams->mollyIter;
  int i, iter;
  HGArrayBigReal dsq;
  BigReal tmp;
  HGArrayInt ial, ibl;
  HGArrayVector ref;  // reference position
  HGArrayVector refab;  // reference vector
  HGArrayBigReal rmass;  // 1 / mass
  BigReal *lambda;  // Lagrange multipliers
  CompAtom *avg;  // averaged position
  int numLambdas = 0;
  HGArrayInt fixed;  // is atom fixed?

  //  iout<<iINFO << "mollyAverage: "<<std::endl<<endi;
  p_avg.resize(numAtoms);
  for ( i=0; i<numAtoms; ++i ) p_avg[i] = p[i];

  for ( int ig = 0; ig < numAtoms; ig += atom[ig].hydrogenGroupSize ) {
    int hgs = atom[ig].hydrogenGroupSize;
    if ( hgs == 1 ) continue;  // only one atom in group
	for ( i = 0; i < hgs; ++i ) {
	  ref[i] = atom[ig+i].position;
	  rmass[i] = 1. / atom[ig+i].mass;
	  fixed[i] = ( simParams->fixedAtomsOn && atom[ig+i].atomFixed );
	  if ( fixed[i] ) rmass[i] = 0.;
	}
	avg = &(p_avg[ig]);
	int icnt = 0;

	if ( ( tmp = mol->rigid_bond_length(atom[ig].id) ) ) {  // for water
	  if ( hgs != 3 ) {
	    NAMD_die("Hydrogen group error caught in mollyAverage().  It's a bug!\n");
	  }
	  if ( !(fixed[1] && fixed[2]) ) {
	    dsq[icnt] = tmp * tmp;  ial[icnt] = 1;  ibl[icnt] = 2;  ++icnt;
	  }
	}
	for ( i = 1; i < hgs; ++i ) {  // normal bonds to mother atom
	  if ( ( tmp = mol->rigid_bond_length(atom[ig+i].id) ) ) {
	    if ( !(fixed[0] && fixed[i]) ) {
	      dsq[icnt] =  tmp * tmp;  ial[icnt] = 0;  ibl[icnt] = i;  ++icnt;
	    }
	  }
	}
	if ( icnt == 0 ) continue;  // no constraints
	numLambdas += icnt;
	molly_lambda.resize(numLambdas);
	lambda = &(molly_lambda[numLambdas - icnt]);
	for ( i = 0; i < icnt; ++i ) {
	  refab[i] = ref[ial[i]] - ref[ibl[i]];
	}
	//	iout<<iINFO<<"hgs="<<hgs<<" m="<<icnt<<std::endl<<endi;
	iter=average(avg,ref,lambda,hgs,icnt,rmass,dsq,ial,ibl,refab,tol,maxiter);
	if ( iter == maxiter ) {
	  iout << iWARN << "Exceeded maximum number of iterations in mollyAverage().\n"<<endi;
	}
  }

  // for ( i=0; i<numAtoms; ++i ) {
  //    if ( ( p_avg[i].position - p[i].position ).length2() > 1.0 ) {
  //      iout << iERROR << "MOLLY moved atom " << (p[i].id + 1) << " from "
  //        << p[i].position << " to " << p_avg[i].position << "\n" << endi;
  //    }
  // }

}


//  MOLLY algorithm part 2
void HomePatch::mollyMollify(Tensor *virial)
{
  Molecule *mol = Node::Object()->molecule;
  SimParameters *simParams = Node::Object()->simParameters;
  Tensor wc;  // constraint virial
  int i;
  HGArrayInt ial, ibl;
  HGArrayVector ref;  // reference position
  CompAtom *avg;  // averaged position
  HGArrayVector refab;  // reference vector
  HGArrayVector force;  // new force
  HGArrayBigReal rmass;  // 1 / mass
  BigReal *lambda;  // Lagrange multipliers
  int numLambdas = 0;
  HGArrayInt fixed;  // is atom fixed?

  for ( int ig = 0; ig < numAtoms; ig += atom[ig].hydrogenGroupSize ) {
    int hgs = atom[ig].hydrogenGroupSize;
    if (hgs == 1 ) continue;  // only one atom in group
	for ( i = 0; i < hgs; ++i ) {
	  ref[i] = atom[ig+i].position;
	  force[i] = f[Results::slow][ig+i];
	  rmass[i] = 1. / atom[ig+i].mass;
	  fixed[i] = ( simParams->fixedAtomsOn && atom[ig+i].atomFixed );
	  if ( fixed[i] ) rmass[i] = 0.;
	}
	int icnt = 0;
	// c-ji I'm only going to mollify water for now
	if ( ( mol->rigid_bond_length(atom[ig].id) ) ) {  // for water
	  if ( hgs != 3 ) {
	    NAMD_die("Hydrogen group error caught in mollyMollify().  It's a bug!\n");
	  }
	  if ( !(fixed[1] && fixed[2]) ) {
	    ial[icnt] = 1;  ibl[icnt] = 2;  ++icnt;
	  }
	}
	for ( i = 1; i < hgs; ++i ) {  // normal bonds to mother atom
	  if ( ( mol->rigid_bond_length(atom[ig+i].id) ) ) {
	    if ( !(fixed[0] && fixed[i]) ) {
	      ial[icnt] = 0;  ibl[icnt] = i;  ++icnt;
	    }
	  }
	}

	if ( icnt == 0 ) continue;  // no constraints
	lambda = &(molly_lambda[numLambdas]);
	numLambdas += icnt;
	for ( i = 0; i < icnt; ++i ) {
	  refab[i] = ref[ial[i]] - ref[ibl[i]];
	}
	avg = &(p_avg[ig]);
	mollify(avg,ref,lambda,force,hgs,icnt,rmass,ial,ibl,refab);
	// store data back to patch
	for ( i = 0; i < hgs; ++i ) {
	  wc += outer(force[i]-f[Results::slow][ig+i],ref[i]);
	  f[Results::slow][ig+i] = force[i];
	}
  }
  // check that there isn't a constant needed here!
  *virial += wc;
  p_avg.resize(0);
}

void HomePatch::checkpoint(void) {
  FullAtomList tmp_a(&atom); checkpoint_atom = tmp_a;
  checkpoint_lattice = lattice;

  // DMK - Atom Separation (water vs. non-water)
  #if NAMD_SeparateWaters != 0
    checkpoint_numWaterAtoms = numWaterAtoms;
  #endif
}

void HomePatch::revert(void) {
  FullAtomList tmp_a(&checkpoint_atom); atom = tmp_a;
  numAtoms = atom.size();
  lattice = checkpoint_lattice;

  // DMK - Atom Separation (water vs. non-water)
  #if NAMD_SeparateWaters != 0
    numWaterAtoms = checkpoint_numWaterAtoms;
  #endif
}

void HomePatch::submitLoadStats(int timestep)
{
  LdbCoordinator::Object()->patchLoad(patchID,numAtoms,timestep);
}


void HomePatch::doPairlistCheck()
{
  SimParameters *simParams = Node::Object()->simParameters;

  if ( ! flags.usePairlists ) {
    flags.pairlistTolerance = 0.;
    flags.maxAtomMovement = 99999.;
    return;
  }

  int i; int n = numAtoms;
  CompAtom *p_i = p.begin();

  if ( flags.savePairlists ) {
    flags.pairlistTolerance = doPairlistCheck_newTolerance;
    flags.maxAtomMovement = 0.;
    doPairlistCheck_newTolerance *= (1. - simParams->pairlistShrink);
    doPairlistCheck_lattice = lattice;
    doPairlistCheck_positions.resize(numAtoms);
    CompAtom *psave_i = doPairlistCheck_positions.begin();
    for ( i=0; i<n; ++i ) { psave_i[i] = p_i[i]; }
    return;
  }

  Lattice &lattice_old = doPairlistCheck_lattice;
  Position center_cur = lattice.unscale(center);
  Position center_old = lattice_old.unscale(center);
  Vector center_delta = center_cur - center_old;
  
  // find max deviation to corner (any neighbor shares a corner)
  BigReal max_cd = 0.;
  for ( i=0; i<2; ++i ) {
    for ( int j=0; j<2; ++j ) {
      for ( int k=0; k<2; ++k ) {
	ScaledPosition corner(	i ? min.x : max.x ,
				j ? min.y : max.y ,
				k ? min.z : max.z );
	Vector corner_delta =
		lattice.unscale(corner) - lattice_old.unscale(corner);
        corner_delta -= center_delta;
	BigReal cd = corner_delta.length2();
        if ( cd > max_cd ) max_cd = cd;
      }
    }
  }
  max_cd = sqrt(max_cd);

  // find max deviation of atoms relative to center
  BigReal max_pd = 0.;
  CompAtom *p_old_i = doPairlistCheck_positions.begin();
  for ( i=0; i<n; ++i ) {
    Vector p_delta = p_i[i].position - p_old_i[i].position;
    p_delta -= center_delta;
    BigReal pd = p_delta.length2();
    if ( pd > max_pd ) max_pd = pd;
  }
  max_pd = sqrt(max_pd);

  BigReal max_tol = max_pd + max_cd;

  flags.maxAtomMovement = max_tol;

  // if ( max_tol > flags.pairlistTolerance ) iout << "tolerance " << max_tol << " > " << flags.pairlistTolerance << "\n" << endi;

  if ( max_tol > ( (1. - simParams->pairlistTrigger) *
				doPairlistCheck_newTolerance ) ) {
    doPairlistCheck_newTolerance *= (1. + simParams->pairlistGrow);
  }

  if ( max_tol > doPairlistCheck_newTolerance ) {
    doPairlistCheck_newTolerance = max_tol / (1. - simParams->pairlistTrigger);
  }

}

void HomePatch::doGroupSizeCheck()
{
  if ( ! flags.doNonbonded ) return;

  SimParameters *simParams = Node::Object()->simParameters;
  BigReal hgcut = 0.5 * simParams->hgroupCutoff;  hgcut *= hgcut;
  BigReal maxrad2 = 0.;

  FullAtomList::iterator p_i = atom.begin();
  FullAtomList::iterator p_e = atom.end();

  while ( p_i != p_e ) {
    int hgs = p_i->hydrogenGroupSize;
    p_i->nonbondedGroupSize = hgs;
    BigReal x = p_i->position.x;
    BigReal y = p_i->position.y;
    BigReal z = p_i->position.z;
    ++p_i;
    int oversize = 0;
    // limit spatial extent
    for ( int i = 1; i < hgs; ++i ) {
      p_i->nonbondedGroupSize = 0;
      BigReal dx = p_i->position.x - x;
      BigReal dy = p_i->position.y - y;
      BigReal dz = p_i->position.z - z;
      BigReal r2 = dx * dx + dy * dy + dz * dz;
      ++p_i;
      if ( r2 > hgcut ) oversize = 1;
      else if ( r2 > maxrad2 ) maxrad2 = r2;
    }
    // also limit to at most 4 atoms per group
    if ( oversize || hgs > 4 ) {
      p_i -= hgs;
      for ( int i = 0; i < hgs; ++i ) {
        p_i->nonbondedGroupSize = 1;
        ++p_i;
      }
    }
  }

  flags.maxGroupRadius = sqrt(maxrad2);

}

void HomePatch::doMarginCheck()
{
  SimParameters *simParams = Node::Object()->simParameters;

  BigReal sysdima = lattice.a_r().unit() * lattice.a();
  BigReal sysdimb = lattice.b_r().unit() * lattice.b();
  BigReal sysdimc = lattice.c_r().unit() * lattice.c();

  BigReal minSize = simParams->patchDimension - simParams->margin;

  if ( ( (max.x - min.x)*aAway*sysdima < minSize*0.9999 ) ||
       ( (max.y - min.y)*bAway*sysdimb < minSize*0.9999 ) ||
       ( (max.z - min.z)*cAway*sysdimc < minSize*0.9999 ) ) {

    NAMD_die("Periodic cell has become too small for original patch grid!\n"
      "Possible solutions are to restart from a recent checkpoint,\n"
      "increase margin, or disable useFlexibleCell for liquid simulation.");
  }

  BigReal cutoff = simParams->cutoff;

  BigReal margina = 0.5 * ( (max.x - min.x) * aAway - cutoff / sysdima );
  BigReal marginb = 0.5 * ( (max.y - min.y) * bAway - cutoff / sysdimb );
  BigReal marginc = 0.5 * ( (max.z - min.z) * cAway - cutoff / sysdimc );

  if ( (margina < -0.0001) || (marginb < -0.0001) || (marginc < -0.0001) ) {
    NAMD_die("Periodic cell has become too small for original patch grid!\n"
      "There are probably many margin violations already reported.\n"
      "Possible solutions are to restart from a recent checkpoint,\n"
      "increase margin, or disable useFlexibleCell for liquid simulation.");
  }

  BigReal minx = min.x - margina;
  BigReal miny = min.y - marginb;
  BigReal minz = min.z - marginc;
  BigReal maxx = max.x + margina;
  BigReal maxy = max.y + marginb;
  BigReal maxz = max.z + marginc;

  int xdev, ydev, zdev;
  int problemCount = 0;

  FullAtomList::iterator p_i = atom.begin();
  FullAtomList::iterator p_e = atom.end();
  for ( ; p_i != p_e; ++p_i ) {

    ScaledPosition s = lattice.scale(p_i->position);

    // check if atom is within bounds
    if (s.x < minx) xdev = 0;
    else if (maxx <= s.x) xdev = 2; 
    else xdev = 1;

    if (s.y < miny) ydev = 0;
    else if (maxy <= s.y) ydev = 2; 
    else ydev = 1;

    if (s.z < minz) zdev = 0;
    else if (maxz <= s.z) zdev = 2; 
    else zdev = 1;

    if (mInfo[xdev][ydev][zdev]) { // somewhere else to be
	++problemCount;
    }

  }

  marginViolations = problemCount;
  // if ( problemCount ) {
  //     iout << iERROR <<
  //       "Found " << problemCount << " margin violations!\n" << endi;
  // } 

}


void
HomePatch::doAtomMigration()
{
  int i;

  for (i=0; i<numNeighbors; i++) {
    realInfo[i].mList.resize(0);
  }

  // Purge the AtomMap
  AtomMap::Object()->unregisterIDs(patchID,pExt.begin(),pExt.end());

  // Determine atoms that need to migrate

  BigReal minx = min.x;
  BigReal miny = min.y;
  BigReal minz = min.z;
  BigReal maxx = max.x;
  BigReal maxy = max.y;
  BigReal maxz = max.z;

  int xdev, ydev, zdev;
  int delnum = 0;

  FullAtomList::iterator atom_i = atom.begin();
  FullAtomList::iterator atom_e = atom.end();

  // DMK - Atom Separation (water vs. non-water)
  #if NAMD_SeparateWaters != 0
    FullAtomList::iterator atom_first = atom_i;
    int numLostWaterAtoms = 0;
  #endif

  while ( atom_i != atom_e ) {
    if ( atom_i->hydrogenGroupSize ) {

      ScaledPosition s = lattice.scale(atom_i->position);

      // check if atom is within bounds
      if (s.x < minx) xdev = 0;
      else if (maxx <= s.x) xdev = 2;
      else xdev = 1;

      if (s.y < miny) ydev = 0;
      else if (maxy <= s.y) ydev = 2;
      else ydev = 1;

      if (s.z < minz) zdev = 0;
      else if (maxz <= s.z) zdev = 2;
      else zdev = 1;

    }

    if (mInfo[xdev][ydev][zdev]) { // process atom for migration
                                    // Don't migrate if destination is myself

      // See if we have a migration list already
      MigrationList &mCur = mInfo[xdev][ydev][zdev]->mList;
      DebugM(3,"Migrating atom " << atom_i->id << " from patch "
		<< patchID << " with position " << atom_i->position << "\n");
      mCur.add(*atom_i);

      ++delnum;


      // DMK - Atom Separation (water vs. non-water)
      #if NAMD_SeparateWaters != 0
        // Check to see if this atom is part of a water molecule.  If
        //   so, numWaterAtoms needs to be adjusted to refect the lost of
        //   this atom.
        // NOTE: The atom separation code assumes that if the oxygen
        //   atom of the hydrogen group making up the water molecule is
        //   migrated to another HomePatch, the hydrogens will also
        //   move!!!
        int atomIndex = atom_i - atom_first;
        if (atomIndex < numWaterAtoms)
          numLostWaterAtoms++;
      #endif


    } else {

      if ( delnum ) { *(atom_i-delnum) = *atom_i; }

    }

    ++atom_i;
  }

  // DMK - Atom Separation (water vs. non-water)
  #if NAMD_SeparateWaters != 0
    numWaterAtoms -= numLostWaterAtoms;
  #endif


  int delpos = numAtoms - delnum;
  DebugM(4,"numAtoms " << numAtoms << " deleted " << delnum << "\n");
  atom.del(delpos,delnum);

  numAtoms = atom.size();

  PatchMgr::Object()->sendMigrationMsgs(patchID, realInfo, numNeighbors);

  // signal depositMigration() that we are inMigration mode
  inMigration = true;

  // Drain the migration message buffer
  for (i=0; i<numMlBuf; i++) {
     DebugM(1, "Draining migration buffer ("<<i<<","<<numMlBuf<<")\n");
     depositMigration(msgbuf[i]);
  }
  numMlBuf = 0;
     
  if (!allMigrationIn) {
    DebugM(3,"All Migrations NOT in, we are suspending patch "<<patchID<<"\n");
    migrationSuspended = true;
    sequencer->suspend();
    migrationSuspended = false;
  }
  allMigrationIn = false;

  inMigration = false;
  marginViolations = 0;
}

void 
HomePatch::depositMigration(MigrateAtomsMsg *msg)
{

  if (!inMigration) { // We have to buffer changes due to migration
		      // until our patch is in migration mode
    msgbuf[numMlBuf++] = msg;
    return;
  } 


  // DMK - Atom Separation (water vs. non-water)
  #if NAMD_SeparateWaters != 0


    // Merge the incoming list of atoms with the current list of
    //   atoms.  Note that mergeSeparatedAtomList() will apply any
    //   required transformations to the incoming atoms as it is
    //   separating them.
    mergeAtomList(msg->migrationList);


  #else

    // Merge the incoming list of atoms with the current list of
    // atoms.  Apply transformations to the incoming atoms as they are
    // added to this patch's list.
    {
      MigrationList &migrationList = msg->migrationList;
      MigrationList::iterator mi;
      Transform mother_transform;
      for (mi = migrationList.begin(); mi != migrationList.end(); mi++) {
        DebugM(1,"Migrating atom " << mi->id << " to patch "
		  << patchID << " with position " << mi->position << "\n"); 
        if ( mi->hydrogenGroupSize ) {
          mi->position = lattice.nearest(mi->position,center,&(mi->transform));
          mother_transform = mi->transform;
        } else {
          mi->position = lattice.reverse_transform(mi->position,mi->transform);
          mi->position = lattice.apply_transform(mi->position,mother_transform);
          mi->transform = mother_transform;
        }
        atom.add(*mi);
      }
    }


  #endif // if (NAMD_SeparateWaters != 0)


  numAtoms = atom.size();
  delete msg;

  DebugM(3,"Counter on " << patchID << " = " << patchMigrationCounter << "\n");
  if (!--patchMigrationCounter) {
    DebugM(3,"All Migrations are in for patch "<<patchID<<"\n");
    allMigrationIn = true;
    patchMigrationCounter = numNeighbors;
    if (migrationSuspended) {
      DebugM(3,"patch "<<patchID<<" is being awakened\n");
      migrationSuspended = false;
      sequencer->awaken();
      return;
    }
    else {
       DebugM(3,"patch "<<patchID<<" was not suspended\n");
    }
  }
}



// DMK - Atom Separation (water vs. non-water)
#if NAMD_SeparateWaters != 0

// This function will separate waters from non-waters in the current
//   atom list (regardless of whether or not the atom list is has been
//   sorted yet or not).
void HomePatch::separateAtoms() {

  // Basic Idea:  Iterate through all the atoms in the current list
  //   of atoms.  Pack the waters in the current atoms list and move
  //   the non-waters to the scratch list.  Once the atoms have all
  //   been separated, copy the non-waters to the end of the waters.
  // NOTE:  This code does not assume that the atoms list has been
  //   separated in any manner.

  // NOTE: Sanity check - Doesn't look like the default constructor actually
  //   adds any atoms but does set numAtoms. ???
  if (atom.size() < 0) return;  // Nothing to do.

  // Resize the scratch FullAtomList (tempAtom)
  tempAtom.resize(numAtoms);  // NOTE: Worst case: all non-water

  // Define size of a water hydrogen group
  int wathgsize = 3;
  if (simParams->watmodel == WAT_TIP4) wathgsize = 4;
  else if (simParams->watmodel == WAT_SWM4) wathgsize = 5;

  // Iterate through all the atoms
  int i = 0;
  int waterIndex = 0;
  int nonWaterIndex = 0;
  while (i < numAtoms) {

    FullAtom &atom_i = atom[i];
    Mass mass = atom_i.mass;
    int hgs = atom_i.hydrogenGroupSize; 
    // Check to see if this hydrogen group is a water molecule
    if (IS_HYDROGEN_GROUP_WATER(hgs, mass)) {

      // Move this hydrogen group up in the current atom list
      if (waterIndex != i) {
        atom[waterIndex    ] = atom[i    ];  // Oxygen
        atom[waterIndex + 1] = atom[i + 1];  // Hydrogen
        atom[waterIndex + 2] = atom[i + 2];  // Hydrogen
        if (wathgsize > 3) atom[waterIndex + 3] = atom[i + 3];  // lonepair
        if (wathgsize > 4) atom[waterIndex + 4] = atom[i + 4];  // drude
          // actual Drude water is arranged:  O D LP H H
      }

      waterIndex += wathgsize;
      i += wathgsize;

    } else {

      // Move this hydrogen group into non-water (scratch) atom list
      for (int j = 0; j < hgs; j++)
        tempAtom[nonWaterIndex + j] = atom[i + j];

      nonWaterIndex += hgs;
      i += hgs;
    }

  } // end iterating through atoms

  // Iterate through the non-water (scratch) atom list, adding the
  //   atoms to the end of the atom list.
  // NOTE: This could be done with a straight memcpy if the internal
  //   data structures of ResizeArray could be accessed directly.
  //   Or, perhaps add a member to ResizeArray that can add a consecutive
  //   list of elements starting at a particular index (would be cleaner).
  for (i = 0; i < nonWaterIndex; i++)
    atom[waterIndex + i] = tempAtom[i];

  // Set numWaterAtoms
  numWaterAtoms = waterIndex;
}


// This function will merge the given list of atoms (not assumed to
//   be separated) with the current list of atoms (already assumed
//   to be separated).
// NOTE: This function applies the transformations to the incoming
//   atoms as it is separating them.
void HomePatch::mergeAtomList(FullAtomList &al) {

  // Sanity check
  if (al.size() <= 0) return;  // Nothing to do

  const int orig_atomSize = atom.size();
  const int orig_alSize = al.size();

  // Resize the atom list (will eventually hold contents of both lists)
  atom.resize(orig_atomSize + orig_alSize); // NOTE: Will have contents of both


  #if 0  // version where non-waters are moved to scratch first

  
  // Basic Idea:  The current list is separated already so copy the
  //   non-water atoms out of it into the scratch atom array.  Then
  //   separate the incoming/given list (al), adding the waters to the
  //   end of the waters in atom list and non-waters to the end of the
  //   scratch list.  At this point, all waters are in atom list and all
  //   non-waters are in the scratch list so just copy the scratch list
  //   to the end of the atom list.
  // NOTE: If al is already separated and the number of waters in it
  //   is know, could simply move the non-waters in atoms back by that
  //   amount and directly copy the waters in al into the created gap
  //   and the non-waters in al to the end.  Leave this as an
  //   optimization for later since I'm not sure if this would actually
  //   do better as the combining code (for combining migration
  //   messages) would also have to merge the contents of the atom lists
  //   they carry.  Generally speaking, there is probably a faster way
  //   to do this, but this will get it working.

  // Copy all the non-waters in the current atom list into the
  //   scratch atom list.
  const int orig_atom_numNonWaters = orig_atomSize - numWaterAtoms;
  tempAtom.resize(orig_atom_numNonWaters + al.size()); // NOTE: Worst case
  for (int i = 0; i < orig_atom_numNonWaters; i++)
    tempAtom[i] = atom[numWaterAtoms + i];

  // Separate the contents of the given atom list (applying the
  // transforms as needed)
  int atom_waterIndex = numWaterAtoms;
  int atom_nonWaterIndex = orig_atom_numNonWaters;
  int i = 0;
  while (i < orig_alSize) {

    FullAtom &atom_i = al[i];
    int hgs = atom_i.hydrogenGroupSize;
    Mass mass = atom_i.mass;

    if (IS_HYDROGEN_GROUP_WATER(hgs, mass)) {

      // Apply the transforms

      // Oxygen (@ +0)
      al[i].position = lattice.nearest(al[i].position, center, &(al[i].transform));
      Transform mother_transform = al[i].transform;

      // Hydrogen (@ +1)
      al[i+1].position = lattice.reverse_transform(al[i+1].position, al[i+1].transform);
      al[i+1].position = lattice.apply_transform(al[i+1].position, mother_transform);
      al[i+1].transform = mother_transform;

      // Hydrogen (@ +2)
      al[i+2].position = lattice.reverse_transform(al[i+2].position, al[i+2].transform);
      al[i+2].position = lattice.apply_transform(al[i+2].position, mother_transform);
      al[i+2].transform = mother_transform;

      // Add to the end of the waters in the current list of atoms
      atom[atom_waterIndex    ] = al[i    ];
      atom[atom_waterIndex + 1] = al[i + 1];
      atom[atom_waterIndex + 2] = al[i + 2];

      atom_waterIndex += 3;
      i += 3;

    } else {

      // Apply the transforms

      // Non-Hydrogen (@ +0)
      al[i].position = lattice.nearest(al[i].position, center, &(al[i].transform));
      Transform mother_transform = al[i].transform;

      // Hydrogens (@ +1 -> +(hgs-1))
      for (int j = 1; j < hgs; j++) {
        al[i+j].position = lattice.reverse_transform(al[i+j].position, al[i+j].transform);
        al[i+j].position = lattice.apply_transform(al[i+j].position, mother_transform);
        al[i+j].transform = mother_transform;
      }

      // Add to the end of the non-waters (scratch) atom list
      for (int j = 0; j < hgs; j++)
        tempAtom[atom_nonWaterIndex + j] = al[i + j];

      atom_nonWaterIndex += hgs;
      i += hgs;
    }

  } // end while iterating through given atom list

  // Copy all the non-waters to the end of the current atom list
  for (int i = 0; i < atom_nonWaterIndex; i++)
    atom[atom_waterIndex + i] = tempAtom[i];

  // Set numWaterAtoms and numAtoms
  numWaterAtoms = atom_waterIndex;
  numAtoms = atom.size();


  #else


  // Basic Idea:  Count the number of water atoms in the incoming atom
  //   list then move the non-waters back in the current atom list to
  //   make room for the incoming waters.  Once there is room in the
  //   current list, separate the incoming list as the atoms are being
  //   added to the current list.
  // NOTE:  Since the incoming atom list is likely to be small,
  //   iterating over its hydrogen groups twice should not be too bad.
  // NOTE:  This code assumes the current list is already separated,
  //   the incoming list may not be separated, and the transforms are
  //   applied to the incoming atoms as the separation occurs.

  // size of a water hydrogen group
  int wathgsize = 3;
  if (simParams->watmodel == WAT_TIP4) wathgsize = 4;
  else if (simParams->watmodel == WAT_SWM4) wathgsize = 5;

  // Count the number of waters in the given atom list
  int al_numWaterAtoms = 0;
  int i = 0;
  while (i < orig_alSize) {

    FullAtom &atom_i = al[i];
    int hgs = atom_i.hydrogenGroupSize;
    Mass mass = atom_i.mass;

    if (IS_HYDROGEN_GROUP_WATER(hgs, mass)) {
      al_numWaterAtoms += wathgsize;
    }

    i += hgs;
  }

  // Move all of the non-waters in the current atom list back (to a
  //   higher index) by the number of waters in the given list.
  if (al_numWaterAtoms > 0) {
    for (i = orig_atomSize - 1; i >= numWaterAtoms; i--) {
      atom[i + al_numWaterAtoms] = atom[i];
    }
  }

  // Separte the atoms in the given atom list.  Perform the
  //   transformations on them and then add them to the appropriate
  //   location in the current atom list.
  int atom_waterIndex = numWaterAtoms;
  int atom_nonWaterIndex = orig_atomSize + al_numWaterAtoms;
  i = 0;
  while (i < orig_alSize) {

    FullAtom &atom_i = al[i];
    int hgs = atom_i.hydrogenGroupSize;
    Mass mass = atom_i.mass;

    if (IS_HYDROGEN_GROUP_WATER(hgs, mass)) {

      // Apply the transforms

      // Oxygen (@ +0)
      al[i].position = lattice.nearest(al[i].position, center, &(al[i].transform));
      Transform mother_transform = al[i].transform;

      // Hydrogen (@ +1)
      al[i+1].position = lattice.reverse_transform(al[i+1].position, al[i+1].transform);
      al[i+1].position = lattice.apply_transform(al[i+1].position, mother_transform);
      al[i+1].transform = mother_transform;

      // Hydrogen (@ +2)
      al[i+2].position = lattice.reverse_transform(al[i+2].position, al[i+2].transform);
      al[i+2].position = lattice.apply_transform(al[i+2].position, mother_transform);
      al[i+2].transform = mother_transform;

      // Add to the end of the waters in the current list of atoms
      atom[atom_waterIndex    ] = al[i    ];
      atom[atom_waterIndex + 1] = al[i + 1];
      atom[atom_waterIndex + 2] = al[i + 2];

      if (wathgsize > 3) atom[atom_waterIndex + 3] = al[i + 3];

      atom_waterIndex += wathgsize;
      i += wathgsize;

    } else {

      // Apply the transforms

      // Non-Hydrogen (@ +0)
      al[i].position = lattice.nearest(al[i].position, center, &(al[i].transform));
      Transform mother_transform = al[i].transform;

      // Hydrogens (@ +1 -> +(hgs-1))
      for (int j = 1; j < hgs; j++) {
        al[i+j].position = lattice.reverse_transform(al[i+j].position, al[i+j].transform);
        al[i+j].position = lattice.apply_transform(al[i+j].position, mother_transform);
        al[i+j].transform = mother_transform;
      }

      // Add to the end of the non-waters (scratch) atom list
      for (int j = 0; j < hgs; j++)
        atom[atom_nonWaterIndex + j] = al[i + j];

      atom_nonWaterIndex += hgs;
      i += hgs;
    }

  } // end while iterating through given atom list

  // Set numWaterAtoms and numAtoms
  numWaterAtoms = atom_waterIndex;
  numAtoms = atom_nonWaterIndex;

  #endif
}

#endif



inline void lubksb(HGMatrixBigReal &a, int n, HGArrayInt &indx,
                                              HGArrayBigReal &b)
{
	int i,ii=-1,ip,j;
	double sum;

	for (i=0;i<n;i++) {
		ip=indx[i];
		sum=b[ip];
		b[ip]=b[i];
		if (ii >= 0)
			for (j=ii;j<i;j++) sum -= a[i][j]*b[j];
		else if (sum) ii=i;
		b[i]=sum;
	}
	for (i=n-1;i>=0;i--) {
		sum=b[i];
		for (j=i+1;j<n;j++) sum -= a[i][j]*b[j];
		b[i]=sum/a[i][i];
	}
}


inline void ludcmp(HGMatrixBigReal &a, int n, HGArrayInt &indx, BigReal *d)
{

	int i,imax,j,k;
	double big,dum,sum,temp;
	HGArrayBigReal vv;
	*d=1.0;
	for (i=0;i<n;i++) {
		big=0.0;
		for (j=0;j<n;j++)
			if ((temp=fabs(a[i][j])) > big) big=temp;
		if (big == 0.0) NAMD_die("Singular matrix in routine ludcmp\n");
		vv[i]=1.0/big;
	}
	for (j=0;j<n;j++) {
		for (i=0;i<j;i++) {
			sum=a[i][j];
			for (k=0;k<i;k++) sum -= a[i][k]*a[k][j];
			a[i][j]=sum;
		}
		big=0.0;
		for (i=j;i<n;i++) {
			sum=a[i][j];
			for (k=0;k<j;k++)
				sum -= a[i][k]*a[k][j];
			a[i][j]=sum;
			if ( (dum=vv[i]*fabs(sum)) >= big) {
				big=dum;
				imax=i;
			}
		}
		if (j != imax) {
			for (k=0;k<n;k++) {
				dum=a[imax][k];
				a[imax][k]=a[j][k];
				a[j][k]=dum;
			}
			*d = -(*d);
			vv[imax]=vv[j];
		}
		indx[j]=imax;
		if (a[j][j] == 0.0) a[j][j]=TINY;
		if (j != n-1) {
			dum=1.0/(a[j][j]);
			for (i=j+1;i<n;i++) a[i][j] *= dum;
		}
	}
}


inline void G_q(const HGArrayVector &refab, HGMatrixVector &gqij,
     const int n, const int m, const HGArrayInt &ial, const HGArrayInt &ibl) {
  int i; 
  // step through the rows of the matrix
  for(i=0;i<m;i++) {
    gqij[i][ial[i]]=2.0*refab[i];
    gqij[i][ibl[i]]=-gqij[i][ial[i]];
  }
};


// c-ji code for MOLLY 7-31-99
int average(CompAtom *qtilde,const HGArrayVector &q,BigReal *lambda,const int n,const int m, const HGArrayBigReal &imass, const HGArrayBigReal &length2, const HGArrayInt &ial, const HGArrayInt &ibl, const HGArrayVector &refab, const BigReal tolf, const int ntrial) {
  //  input:  n = length of hydrogen group to be averaged (shaked)
  //          q[n] = vector array of original positions
  //          m = number of constraints
  //          imass[n] = inverse mass for each atom
  //          length2[m] = square of reference bond length for each constraint
  //          ial[m] = atom a in each constraint 
  //          ibl[m] = atom b in each constraint 
  //          refab[m] = vector of q_ial(i) - q_ibl(i) for each constraint
  //          tolf = function error tolerance for Newton's iteration
  //          ntrial = max number of Newton's iterations
  //  output: lambda[m] = double array of lagrange multipliers (used by mollify)
  //          qtilde[n] = vector array of averaged (shaked) positions

  int k,k1,i,j;
  BigReal errx,errf,d,tolx;

  HGArrayInt indx;
  HGArrayBigReal p;
  HGArrayBigReal fvec;
  HGMatrixBigReal fjac;
  HGArrayVector avgab;
  HGMatrixVector grhs;
  HGMatrixVector auxrhs;
  HGMatrixVector glhs;

  //  iout <<iINFO << "average: n="<<n<<" m="<<m<<std::endl<<endi;
  tolx=tolf; 
  
  // initialize lambda, globalGrhs

  for (i=0;i<m;i++) {
    lambda[i]=0.0;
  }

  // define grhs, auxrhs for all iterations
  // grhs= g_x(q)
  //
  G_q(refab,grhs,n,m,ial,ibl);
  for (k=1;k<=ntrial;k++) {
    //    usrfun(qtilde,q0,lambda,fvec,fjac,n,water); 
    HGArrayBigReal gij;
    // this used to be main loop of usrfun
    // compute qtilde given q0, lambda, IMASSes
    {
      BigReal multiplier;
      HGArrayVector tmp;
      for (i=0;i<m;i++) {
	multiplier = lambda[i];
	// auxrhs = M^{-1}grhs^{T}
	for (j=0;j<n;j++) {
	  auxrhs[i][j]=multiplier*imass[j]*grhs[i][j];
	}
      }
      for (j=0;j<n;j++) {
	//      tmp[j]=0.0;      
	for (i=0;i<m;i++) {
	  tmp[j]+=auxrhs[i][j];
	}
      }
 
      for (j=0;j<n;j++) {
	qtilde[j].position=q[j]+tmp[j];
      }
      //      delete [] tmp;
    }
  
    for ( i = 0; i < m; i++ ) {
      avgab[i] = qtilde[ial[i]].position - qtilde[ibl[i]].position;
    }

    //  iout<<iINFO << "Calling Jac" << std::endl<<endi;
    //  Jac(qtilde, q0, fjac,n,water);
    {
      //  Vector glhs[3*n+3];

      HGMatrixVector grhs2;

      G_q(avgab,glhs,n,m,ial,ibl);
#ifdef DEBUG0
      iout<<iINFO << "G_q:" << std::endl<<endi;
      for (i=0;i<m;i++) {
	iout<<iINFO << glhs[i*n+0] << " " << glhs[i*n+1] << " " << glhs[i*n+2] << std::endl<<endi;
      }
#endif
      //      G_q(refab,grhs2,m,ial,ibl);
      // update with the masses
      for (j=0; j<n; j++) { // number of atoms
	for (i=0; i<m;i++) { // number of constraints
	  grhs2[i][j] = grhs[i][j]*imass[j];
	}
      }

      // G_q(qtilde) * M^-1 G_q'(q0) =
      // G_q(qtilde) * grhs'
      for (i=0;i<m;i++) { // number of constraints
	for (j=0;j<m;j++) { // number of constraints
	  fjac[i][j] = 0; 
	  for (k1=0;k1<n;k1++) {
	    fjac[i][j] += glhs[i][k1]*grhs2[j][k1]; 
	  }
	}
      }
#ifdef DEBUG0  
      iout<<iINFO << "glhs" <<endi;
      for(i=0;i<9;i++) {
	iout<<iINFO << glhs[i] << ","<<endi;
      }
      iout<<iINFO << std::endl<<endi;
      for(i=0;i<9;i++) {
	iout<<iINFO << grhs2[i] << ","<<endi;
      }
      iout<<iINFO << std::endl<<endi;
#endif
      //      delete[] grhs2;
    }
    // end of Jac calculation
#ifdef DEBUG0
    iout<<iINFO << "Jac" << std::endl<<endi;
    for (i=0;i<m;i++) 
      for (j=0;j<m;j++)
	iout<<iINFO << fjac[i][j] << " "<<endi;
    iout<< std::endl<<endi;
#endif
    // calculate constraints in gij for n constraints this being a water
    //  G(qtilde, gij, n, water);
    for (i=0;i<m;i++) {
      gij[i]=avgab[i]*avgab[i]-length2[i];
    }
#ifdef DEBUG0
    iout<<iINFO << "G" << std::endl<<endi;
    iout<<iINFO << "( "<<endi;
    for(i=0;i<m-1;i++) {
      iout<<iINFO << gij[i] << ", "<<endi;
    }
    iout<<iINFO << gij[m-1] << ")" << std::endl<<endi;
#endif
    // fill the return vector
    for(i=0;i<m;i++) {
      fvec[i] = gij[i];
    }
    // free up the constraints
    //    delete[] gij;
    // continue Newton's iteration    
    errf=0.0;
    for (i=0;i<m;i++) errf += fabs(fvec[i]);
#ifdef DEBUG0
    iout<<iINFO << "errf: " << errf << std::endl<<endi;
#endif
    if (errf <= tolf) {
      break;
    }
    for (i=0;i<m;i++) p[i] = -fvec[i];
    //    iout<<iINFO << "Doing dcmp in average " << std::endl<<endi;
    ludcmp(fjac,m,indx,&d);
    lubksb(fjac,m,indx,p);

    errx=0.0;
    for (i=0;i<m;i++) {
      errx += fabs(p[i]);
    }
    for (i=0;i<m;i++)  
      lambda[i] += p[i];

#ifdef DEBUG0
    iout<<iINFO << "lambda:" << lambda[0] 
	 << " " << lambda[1] << " " << lambda[2] << std::endl<<endi;
    iout<<iINFO << "errx: " << errx << std::endl<<endi;
#endif
    if (errx <= tolx) break;
#ifdef DEBUG0
    iout<<iINFO << "Qtilde:" << std::endl<<endi;
    iout<<iINFO << qtilde[0].position << " " << qtilde[1].position << " " << qtilde[2].position << std::endl<<endi; 
#endif
  }
#ifdef DEBUG
  iout<<iINFO << "LAMBDA:" << lambda[0] << " " << lambda[1] << " " << lambda[2] << std::endl<<endi;
#endif

  return k; // 
}

void mollify(CompAtom *qtilde,const HGArrayVector &q0,const BigReal *lambda, HGArrayVector &force,const int n, const int m, const HGArrayBigReal &imass,const HGArrayInt &ial,const HGArrayInt &ibl,const HGArrayVector &refab) {
  int i,j,k;
  BigReal d;
  HGMatrixBigReal fjac;
  Vector zero(0.0,0.0,0.0);
  
  HGArrayVector tmpforce;
  HGArrayVector tmpforce2;
  HGArrayVector y;
  HGMatrixVector grhs;
  HGMatrixVector glhs;
  HGArrayBigReal aux;
  HGArrayInt indx;

  for(i=0;i<n;i++) {
    tmpforce[i]=imass[i]*force[i];
  }

  HGMatrixVector grhs2;
  HGArrayVector avgab;

  for ( i = 0; i < m; i++ ) {
	avgab[i] = qtilde[ial[i]].position - qtilde[ibl[i]].position;
  }

  G_q(avgab,glhs,n,m,ial,ibl);
  G_q(refab,grhs,n,m,ial,ibl);
  // update with the masses
  for (j=0; j<n; j++) { // number of atoms
	for (i=0; i<m;i++) { // number of constraints
	  grhs2[i][j] = grhs[i][j]*imass[j];
	}
  }

  // G_q(qtilde) * M^-1 G_q'(q0) =
  // G_q(qtilde) * grhs'
  for (i=0;i<m;i++) { // number of constraints
	for (j=0;j<m;j++) { // number of constraints
	  fjac[j][i] = 0; 
	  for (k=0;k<n;k++) {
	    fjac[j][i] += glhs[i][k]*grhs2[j][k]; 
	  }
	}
  }

  // aux=gqij*tmpforce
  //  globalGrhs::computeGlobalGrhs(q0,n,water);
  //  G_q(refab,grhs,m,ial,ibl);
  for(i=0;i<m;i++) {
    aux[i]=0.0;
    for(j=0;j<n;j++) {
      aux[i]+=grhs[i][j]*tmpforce[j];
    }
  }

  ludcmp(fjac,m,indx,&d);
  lubksb(fjac,m,indx,aux);

  for(j=0;j<n;j++) {
    y[j] = zero;
    for(i=0;i<m;i++) {
      y[j] += aux[i]*glhs[i][j];
    }
  }
  for(i=0;i<n;i++) {
    y[i]=force[i]-y[i];
  }
    
  // gqq12*y
  for(i=0;i<n;i++) {
    tmpforce2[i]=imass[i]*y[i];
  }

  // here we assume that tmpforce is initialized to zero.
  for (i=0;i<n;i++) {
    tmpforce[i]=zero;
  }
  
  for (j=0;j<m;j++) {
    Vector tmpf = 2.0*lambda[j]*(tmpforce2[ial[j]]-tmpforce2[ibl[j]]);
    tmpforce[ial[j]] += tmpf;
    tmpforce[ibl[j]] -= tmpf;
  }
  // c-ji the other bug for 2 constraint water was this line (2-4-99)
  //  for(i=0;i<m;i++) {
  for(i=0;i<n;i++) {
    force[i]=tmpforce[i]+y[i];
  }

}

#if CMK_PERSISTENT_COMM
void HomePatch::destoryPersistComm()
{
     for (int i=0; i<nphs; i++) {
       CmiDestoryPersistent(localphs[i]);
     }
     phsReady = 0;
}
#endif
