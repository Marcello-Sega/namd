/**
***  Copyright (c) 1995, 1996, 1997, 1998, 1999, 2000 by
***  The Board of Trustees of the University of Illinois.
***  All rights reserved.
**/

/*
   Toplevel routines for initializing a Node for a simulation
   one Node per Pe (processor element).
*/

#if !defined(WIN32) || defined(__CYGWIN__)
#include <unistd.h>
#endif
#include "InfoStream.h"
#include "Node.decl.h"
#include "Node.h"
#ifdef DPMTA
#include <pvm3.h>
#endif

#include "ProcessorPrivate.h"

#define MIN_DEBUG_LEVEL 3
//#define DEBUGM
#include "Debug.h"

#include <stdio.h>
#include <converse.h>
#include "memusage.h"
#include "IMDOutput.h"

#ifdef USE_COMM_LIB
#include "ComlibManager.h"
#endif

#include "Lattice.h"
#include "main.decl.h"
#include "main.h"
#include "WorkDistrib.h"
#include "PatchMgr.h"
#include "Patch.h"
#include "Compute.h"
#include "ComputeMap.h"
#include "ComputeMgr.h"
#include "Molecule.h"
#include "HomePatchList.h"
#include "AtomMap.h"
#include "Sequencer.h"
#include "Controller.h"
#include "NamdState.h"
#include "Output.h"
#include "ProxyMgr.h"
#include "PatchMap.h"
#include "PatchMap.inl"
#include "Parameters.h"
#include "SimParameters.h"
#include "Communicate.h"
#include "LdbCoordinator.h"
#include "ScriptTcl.h"
#include "ComputeMgr.decl.h"
#include "ComputePmeMgr.decl.h"
#include "OptPmeMgr.decl.h"
#include "Sync.h"
#include "BackEnd.h"
#include "PDB.h"

#if(CMK_CCS_AVAILABLE && CMK_WEB_MODE)
extern "C" void CApplicationInit();
#endif

#include "DumpBench.h"

#ifdef MEM_OPT_VERSION
#include "CollectionMgr.h"
#include "CollectionMaster.h"
#include "CollectionMgr.decl.h"
#include "CollectionMaster.decl.h"
#endif

#if USE_HPM
extern "C" void HPM_Init(int);
extern "C" void HPM_Start(char *label, int);
extern "C" void HPM_Stop(char *label, int);
extern "C" void HPM_Print(int, int);
#endif

//======================================================================
// Public Functions

//----------------------------------------------------------------------

int eventEndOfTimeStep;
double startupTime;

//----------------------------------------------------------------------
// BOC constructor
Node::Node(GroupInitMsg *msg)
{    
  DebugM(4,"Creating Node\n");
#if(CMK_CCS_AVAILABLE && CMK_WEB_MODE)
  CApplicationInit();
#endif
  if (CkpvAccess(Node_instance) == 0) {
    CkpvAccess(Node_instance) = this;
    eventEndOfTimeStep = traceRegisterUserEvent("EndOfTimeStep");
  } else {
    NAMD_bug("Node::Node() - another instance of Node exists!");
  }

  CkpvAccess(BOCclass_group) = msg->group;
  delete msg;

  CkpvAccess(BOCclass_group).node = thisgroup;

  startupPhase = 0;

  molecule = NULL;
  parameters = NULL;
  simParameters = NULL;
  configList = NULL;
  pdb = NULL;
  state = NULL;
  output = NULL;
  imd = new IMDOutput;

#if USE_HPM
  // assumes that this will be done only on BG/P
  TopoManager *tmgr = new TopoManager();
  int x, y, z;
  tmgr->rankToCoordinates(CkMyPe(), x, y, z, localRankOnNode);
  delete tmgr;
#endif

  DebugM(4,"Creating PatchMap, AtomMap, ComputeMap\n");
  patchMap = PatchMap::Instance();
  atomMap = AtomMap::Instance();
  computeMap = ComputeMap::Instance();

  DebugM(4,"Binding to BOC's\n");
  CProxy_PatchMgr pm(CkpvAccess(BOCclass_group).patchMgr);
  patchMgr = pm.ckLocalBranch();
  CProxy_ProxyMgr prm(CkpvAccess(BOCclass_group).proxyMgr);
  proxyMgr = prm.ckLocalBranch();
  CProxy_WorkDistrib wd(CkpvAccess(BOCclass_group).workDistrib);
  workDistrib = wd.ckLocalBranch();
  CProxy_ComputeMgr cm(CkpvAccess(BOCclass_group).computeMgr);
  computeMgr = cm.ckLocalBranch();
  CProxy_LdbCoordinator lc(CkpvAccess(BOCclass_group).ldbCoordinator);
  ldbCoordinator = lc.ckLocalBranch();

}

//----------------------------------------------------------------------
// ~Node(void) needs to clean up everything.

Node::~Node(void)
{
  delete output;
  delete computeMap;
  delete atomMap;
  delete patchMap;
  delete CkpvAccess(comm);
}

//----------------------------------------------------------------------
// Startup Sequence

void Node::messageStartUp() {
  (CProxy_Node(CkpvAccess(BOCclass_group).node)).startup();
}

void Node::startUp(CkQdMsg *qmsg) {
  delete qmsg;
  (CProxy_Node(CkpvAccess(BOCclass_group).node)).startup();
}

SimParameters *node_simParameters;
Parameters *node_parameters;
Molecule *node_molecule;

extern void registerUserEventsForAllComputeObjs(void);

void Node::startup() {
  int gotoRun = false;
  double newTime;

  if (!CkMyPe()) {
    if (!startupPhase) {
      iout << iINFO << "\n";
      startupTime = CmiWallTimer();
      iout << iINFO << "Entering startup at " << startupTime << " s, ";
    } else {
      newTime = CmiWallTimer();
      iout << iINFO << "Startup phase " << startupPhase-1 << " took "
	   << newTime - startupTime << " s, ";
      startupTime = newTime;
    }
    iout << memusage_MB() << " MB of memory in use\n" << endi;
    fflush(stdout);
  }
  
  switch (startupPhase) {

  case 0:
    #ifdef CHARMIZE_NAMD
    populateAtomDisArrs(startupPhase);
    #endif

    namdOneCommInit(); // Namd1.X style
  break;

  case 1:
    // send & receive molecule, simparameters... (Namd1.X style)
    if (CkMyPe()) {
      namdOneRecv();
    } else {
      namdOneSend();
    }

    #ifdef CHARMIZE_NAMD
    //send charm array proxies
    if(!CkMyPe()){
        AllCharmArrsMsg *arrsMsg = new AllCharmArrsMsg;
        arrsMsg->atomsDis = atomDisArr;
        ((CProxy_Node)thisgroup).sendCharmArrProxies(arrsMsg);        
    }    
    #endif
  break;

  case 2:
    // fix up one-per-node objects
    simParameters = node_simParameters;
    parameters = node_parameters;
    molecule = node_molecule;

    #if USE_HPM
    HPM_Init(localRankOnNode);
    #endif

    #ifdef MEM_OPT_VERSION
    //Allocate CollectionMaster which handles I/O depending on the shiftIOToOne parameter
    if(!CkMyPe()){
	CkChareID collectionMaster;
	if(CkNumPes()>1 && simParameters->shiftIOToOne)
	    collectionMaster = CProxy_CollectionMaster::ckNew(1);
	else
	    collectionMaster = CProxy_CollectionMaster::ckNew(0);
	
	//set CollectionMgr and CollectionMasterHandler's field for CollectionMaster
	CollectionMasterHandler::Object()->setRealMaster(collectionMaster);
	CProxy_CollectionMgr cmgr(CkpvAccess(BOCclass_group).collectionMgr);
	SlaveInitMsg *bcmaster = new SlaveInitMsg;
	bcmaster->master = collectionMaster;
	cmgr.setCollectionMaster(bcmaster);
    }
    #endif

    // take care of inital thread setting
    threadInit();

    // create blank AtomMap
    AtomMap::Object()->allocateMap(molecule->numAtoms);

    if (!CkMyPe()) {
      if (simParameters->useOptPME)
	CkpvAccess(BOCclass_group).computePmeMgr = CProxy_OptPmeMgr::ckNew();
      else 
	CkpvAccess(BOCclass_group).computePmeMgr = CProxy_ComputePmeMgr::ckNew();
    }

  break;

  case 3:     
    if(simParameters->isSendProxySTEnabled()) {
        ProxyMgr::Object()->setSendSpanning();
    }
    if(simParameters->isRecvProxySTEnabled()) {
        ProxyMgr::Object()->setRecvSpanning();
    }
    #ifdef PROCTRACE_DEBUG
    DebugFileTrace::Instance("procTrace");
    #endif

    if (!CkMyPe()) {
      output = new Output; // create output object just on PE(0)
      workDistrib->patchMapInit(); // create space division

      #ifdef MEM_OPT_VERSION
      //create patches without populating them with atoms
      workDistrib->preCreateHomePatches();
      #else      
      workDistrib->createHomePatches(); // load atoms into HomePatch(es)
      #endif
      
      workDistrib->assignNodeToPatch();
      workDistrib->mapComputes();
      ComputeMap::Object()->printComputeMap();

      registerUserEventsForAllComputeObjs();

      workDistrib->sendMaps();
      #ifdef USE_NODEPATCHMGR
      CProxy_NodeProxyMgr npm(CkpvAccess(BOCclass_group).nodeProxyMgr);
      //a node broadcast
      npm.createProxyInfo(PatchMap::Object()->numPatches());
      #endif
    }
    {
        #if defined(NODEAWARE_PROXY_SPANNINGTREE) && defined(USE_NODEPATCHMGR)
        CProxy_NodeProxyMgr npm(CkpvAccess(BOCclass_group).nodeProxyMgr);
        if(CkMyRank()==0) {
            //just need to register once
            npm[CkMyNode()].ckLocalBranch()->registerLocalProxyMgr(CkpvAccess(BOCclass_group).proxyMgr);
        }
        npm[CkMyNode()].ckLocalBranch()->registerLocalPatchMap(CkMyRank(), PatchMap::Object());
        #endif
    }
  break;

  case 4:
    if ( simParameters->PMEOn ) {
      if ( simParameters->useOptPME ) {
	CProxy_OptPmeMgr pme(CkpvAccess(BOCclass_group).computePmeMgr);
	pme[CkMyPe()].initialize(new CkQdMsg);
      }
      else {
	CProxy_ComputePmeMgr pme(CkpvAccess(BOCclass_group).computePmeMgr);
	pme[CkMyPe()].initialize(new CkQdMsg);
      }
    }
    break;
    
  case 5:
    if ( simParameters->PMEOn ) {
      if ( simParameters->useOptPME ) {
	CProxy_OptPmeMgr pme(CkpvAccess(BOCclass_group).computePmeMgr);
	pme[CkMyPe()].initialize_pencils(new CkQdMsg);
      }
      else {
	CProxy_ComputePmeMgr pme(CkpvAccess(BOCclass_group).computePmeMgr);
	pme[CkMyPe()].initialize_pencils(new CkQdMsg);
      }
    }
    if (!CkMyPe()) {
    #ifdef MEM_OPT_VERSION
      workDistrib->initAndSendHomePatch();
    #else
      workDistrib->distributeHomePatches();      
    #endif
    }
  break;

  case 6:
    if ( simParameters->PMEOn ) {
      if ( simParameters->useOptPME ) {
	CProxy_OptPmeMgr pme(CkpvAccess(BOCclass_group).computePmeMgr);
	pme[CkMyPe()].activate_pencils(new CkQdMsg);
      }
      else {
	CProxy_ComputePmeMgr pme(CkpvAccess(BOCclass_group).computePmeMgr);
	pme[CkMyPe()].activate_pencils(new CkQdMsg);
      }
    }
    proxyMgr->createProxies();  // need Home patches before this
    if (!CkMyPe()) LdbCoordinator::Object()->createLoadBalancer();
  break;

  case 7:
    if (!CkMyPe()) {
      ComputeMap::Object()->printComputeMap();
    }
    Sync::Object()->openSync();  // decide if to open local Sync 
    if (proxySendSpanning || proxyRecvSpanning )
      proxyMgr->buildProxySpanningTree();
    DebugM(4,"Creating Computes\n");
    computeMgr->createComputes(ComputeMap::Object());
    DebugM(4,"Building Sequencers\n");
    buildSequencers();
    DebugM(4,"Initializing LDB\n");
    LdbCoordinator::Object()->initialize(patchMap,computeMap);
  break;

  case 8:
    {
	//For debugging
	/*if(!CkMyPe()){
	FILE *dumpFile = fopen("/tmp/NAMD_Bench.dump", "w");
	dumpbench(dumpFile);
	NAMD_die("Normal execution\n");
	}*/
    }
    #ifdef MEM_OPT_VERSION
    if(!CkMyPe()){
	molecule->delEachAtomSigs();
	molecule->delChargeSpace();
	if(!simParameters->freeEnergyOn)
	    molecule->delMassSpace();
	molecule->delOtherEachAtomStructs();
	//we can free the pdb data here to save memory because output
	//is always in the binary format thus saving 3 doubles * numAtoms bytes.
	pdb->delPDBCoreData();
    }
    //decide whether to free memory space for cluster information
    //the condition could be referred to comment for function
    //wrap_coor_int in Output.C
    if(simParameters->wrapAll || simParameters->wrapWater){
	int peOfCollectionMaster = 0;
	if(CkNumPes()>1 && simParameters->shiftIOToOne) peOfCollectionMaster = 1;
	if(CkNumPes()>1 && CkMyPe()!=peOfCollectionMaster)
	    molecule->delClusterSigs();	
    }else{
	molecule->delClusterSigs();
    }
    #endif
    gotoRun = true;
  break;

  default:
    NAMD_bug("Startup Phase has a bug - check case statement");
  break;

  }

  startupPhase++;
  if (!CkMyPe()) {
    if (!gotoRun) {
      CkStartQD(CkIndex_Node::startUp((CkQdMsg*)0),&thishandle);
    } else {
      Node::messageRun();
    }
  }
}

void Node::namdOneCommInit()
{
  if (CkpvAccess(comm) == NULL) {
    CkpvAccess(comm) = new Communicate();
#ifdef DPMTA
    pvmc_init();
#endif
  }
}

#ifdef CHARMIZE_NAMD
//It is only executed on pe(0) and called at startup phase 0
void Node::populateAtomDisArrs(int startupPhase){
    if(CkMyPe()) return;
    if(startupPhase) return;

    int disArrSize = molecule->numAtoms/AtomsDisInfo::ATOMDISNUM;
    int remainAtoms = molecule->numAtoms%AtomsDisInfo::ATOMDISNUM;

    int totalArrSize = disArrSize + (remainAtoms!=0);
    atomDisArr = CProxy_AtomsDisInfo::ckNew(totalArrSize);

    int atomIndex = 0;
    Atom *allAtoms = molecule->getAllAtoms();
    for(int i=0; i<totalArrSize; i++){
        int actualAtomCnt = AtomsDisInfo::ATOMDISNUM;
        if(i==disArrSize) actualAtomCnt = remainAtoms;

        AtomStaticInfoMsg *staticMsg = new(actualAtomCnt, 0)AtomStaticInfoMsg;
        staticMsg->actualNumAtoms = actualAtomCnt;
        for(int j=0; j<actualAtomCnt; j++, atomIndex++){
            staticMsg->atoms[j] = allAtoms[atomIndex];
        }       
        atomDisArr[i].recvStaticInfo(staticMsg);
    }    
}
#endif

// Namd 1.X style Send/Recv of simulation information

void Node::namdOneRecv() {
  if ( CmiMyRank() ) return;

  MIStream *conv_msg;

  // Receive molecule and simulation parameter information
  simParameters = node_simParameters = new SimParameters;
  //****** BEGIN CHARMM/XPLOR type changes
  parameters = node_parameters = new Parameters();
  //****** END CHARMM/XPLOR type changes
  molecule = node_molecule = new Molecule(simParameters,parameters);

  DebugM(4, "Getting SimParameters\n");
  conv_msg = CkpvAccess(comm)->newInputStream(0, SIMPARAMSTAG);
  simParameters->receive_SimParameters(conv_msg);

  DebugM(4, "Getting Parameters\n");
  conv_msg = CkpvAccess(comm)->newInputStream(0, STATICPARAMSTAG);
  parameters->receive_Parameters(conv_msg);

  DebugM(4, "Getting Molecule\n");
  conv_msg = CkpvAccess(comm)->newInputStream(0, MOLECULETAG);
  molecule->receive_Molecule(conv_msg);

  DebugM(4, "Done Receiving\n");
}

void Node::namdOneSend() {
  node_simParameters = simParameters;
  node_parameters = parameters;
  node_molecule = molecule;

  MOStream *conv_msg;
  // I'm Pe(0) so I send what I know
  DebugM(4, "Sending SimParameters\n");  
  conv_msg = CkpvAccess(comm)->newOutputStream(ALLBUTME, SIMPARAMSTAG, BUFSIZE);
  simParameters->send_SimParameters(conv_msg);

  DebugM(4, "Sending Parameters\n");
  conv_msg = CkpvAccess(comm)->newOutputStream(ALLBUTME, STATICPARAMSTAG, BUFSIZE);
  parameters->send_Parameters(conv_msg);

  DebugM(4, "Sending Molecule\n");
  int bufSize = BUFSIZE;
  if(molecule->numAtoms>=1000000) bufSize = 16*BUFSIZE;
  conv_msg = CkpvAccess(comm)->newOutputStream(ALLBUTME, MOLECULETAG, bufSize);
  molecule->send_Molecule(conv_msg);
}

void Node::sendCharmArrProxies(AllCharmArrsMsg *msg){
#ifdef CHARMIZE_NAMD
    if(CkMyPe()){
        atomDisArr = msg->atomsDis;        
    }
    delete msg;
#else
    NAMD_die("sendCharmArrProxies should not be called in this case!");
#endif
}

// Initial thread setup

void Node::threadInit() {
  // Thread initialization
  if (CthImplemented()) {
    CthSetStrategyDefault(CthSelf());
  } else {
    NAMD_bug("Node::startup() Oh no, tiny elvis, threads not implemented");
  }
}

//
void Node::buildSequencers() {
  HomePatchList *hpl = PatchMap::Object()->homePatchList();
  ResizeArrayIter<HomePatchElem> ai(*hpl);

  // Controller object is only on Pe(0)
  if ( ! CkMyPe() ) {
    Controller *controller = new Controller(state);
    state->useController(controller);
  }

  // Assign Sequencer to all HomePatch(es)
  for (ai=ai.begin(); ai != ai.end(); ai++) {
    HomePatch *patch = (*ai).patch;
    Sequencer *sequencer = new Sequencer(patch);
    patch->useSequencer(sequencer);
  }
}



//-----------------------------------------------------------------------
// Node run() - broadcast to all nodes
//-----------------------------------------------------------------------
void Node::messageRun() {
  (CProxy_Node(CkpvAccess(BOCclass_group).node)).run();
}


//-----------------------------------------------------------------------
// run(void) runs the specified simulation for the specified number of
// steps, overriding the contents of the configuration file
//-----------------------------------------------------------------------
void Node::run()
{
  // Start Controller (aka scalar Sequencer) on Pe(0)
  if ( ! CkMyPe() ) {
    state->runController();
  }

  DebugM(4, "Starting Sequencers\n");
  // Run Sequencer on each HomePatch - i.e. start simulation
  HomePatchList *hpl = PatchMap::Object()->homePatchList();
  ResizeArrayIter<HomePatchElem> ai(*hpl);
  for (ai=ai.begin(); ai != ai.end(); ai++) {
    HomePatch *patch = (*ai).patch;
    patch->runSequencer();
  }

  if (!CkMyPe()) {
    double newTime = CmiWallTimer();
    iout << iINFO << "Startup phase " << startupPhase-1 << " took "
	 << newTime - startupTime << " s, "
	 << memusage_MB() << " MB of memory in use\n";
    iout << iINFO << "Finished startup at " << newTime << " s, "
	 << memusage_MB() << " MB of memory in use\n\n" << endi;
    fflush(stdout);
  }
  
}


//-----------------------------------------------------------------------
// Node scriptBarrier() - twiddle parameters with simulation halted
//-----------------------------------------------------------------------

void Node::enableScriptBarrier() {
  CkStartQD(CkIndex_Node::scriptBarrier((CkQdMsg*)0),&thishandle);
}

void Node::scriptBarrier(CkQdMsg *qmsg) {
  delete qmsg;
  //script->awaken();
}

void Node::scriptParam(ScriptParamMsg *msg) {
  simParameters->scriptSet(msg->param,msg->value);
  delete msg;
}

void Node::reloadCharges(const char *filename) {
  FILE *file = fopen(filename,"r");
  if ( ! file ) NAMD_die("node::reloadCharges():Error opening charge file.");

  int n = molecule->numAtoms;
  float *charge = new float[n];

  for ( int i = 0; i < n; ++i ) {
    if ( ! fscanf(file,"%f",&charge[i]) )
      NAMD_die("Node::reloadCharges():Not enough numbers in charge file.");
  }

  fclose(file);
  CProxy_Node(thisgroup).reloadCharges(charge,n);
  delete [] charge;
}

void Node::reloadCharges(float charge[], int n) {
  molecule->reloadCharges(charge,n);
}


void Node::sendEnableExitScheduler(void) {
  //CmiPrintf("sendEnableExitScheduler\n");
  CkQdMsg *msg = new CkQdMsg;
  CProxy_Node nodeProxy(thisgroup);
  nodeProxy[0].recvEnableExitScheduler(msg);
}

void Node::recvEnableExitScheduler(CkQdMsg *msg) {
  //CmiPrintf("recvEnableExitScheduler\n");
  delete msg;
  enableExitScheduler();
}

void Node::enableExitScheduler(void) {
  if ( CkMyPe() ) {
    sendEnableExitScheduler();
  } else {
    CkStartQD(CkIndex_Node::exitScheduler((CkQdMsg*)0),&thishandle);
  }
}

void Node::exitScheduler(CkQdMsg *msg) {
  //CmiPrintf("exitScheduler %d\n",CkMyPe());
  CsdExitScheduler();
  delete msg;
}

void Node::sendEnableEarlyExit(void) {
  CkQdMsg *msg = new CkQdMsg;
  CProxy_Node nodeProxy(thisgroup);
  nodeProxy[0].recvEnableEarlyExit(msg);
}

void Node::recvEnableEarlyExit(CkQdMsg *msg) {
  delete msg;
  enableEarlyExit();
}

void Node::enableEarlyExit(void) {
  if ( CkMyPe() ) {
    sendEnableEarlyExit();
  } else {
    CkStartQD(CkIndex_Node::earlyExit((CkQdMsg*)0),&thishandle);
  }
}

void Node::earlyExit(CkQdMsg *msg) {
  iout << iERROR << "Exiting prematurely; see error messages above.\n" << endi;
  BackEnd::exit();
  delete msg;
}


//------------------------------------------------------------------------
// Some odd utilities
//------------------------------------------------------------------------
void Node::saveMolDataPointers(NamdState *state)
{
  this->molecule = state->molecule;
  this->parameters = state->parameters;
  this->simParameters = state->simParameters;
  this->configList = state->configList;
  this->pdb = state->pdb;
  this->state = state;
}

// entry methods for BG/P HPM (performance counters) library
void Node::startHPM() {
#if USE_HPM
  HPM_Start("500 steps", localRankOnNode);
#endif
}

void Node::stopHPM() {
#if USE_HPM
  HPM_Stop("500 steps", localRankOnNode);
  HPM_Print(CkMyPe(), localRankOnNode);
#endif
}

//======================================================================
// Private functions


#include "Node.def.h"

