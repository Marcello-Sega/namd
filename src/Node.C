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
#include <charm++.h>
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
#include "Sync.h"
#include "BackEnd.h"

#if(CMK_CCS_AVAILABLE && CMK_WEB_MODE)
extern "C" void CApplicationInit();
#endif

//======================================================================
// Public Functions

//----------------------------------------------------------------------

int eventEndOfTimeStep;

//----------------------------------------------------------------------
// BOC constructor
Node::Node(GroupInitMsg *msg)
{
  DebugM(4,"Creating Node\n");
#if(CMK_CCS_AVAILABLE && CMK_WEB_MODE)
  CApplicationInit();
#endif
  if (CpvAccess(Node_instance) == 0) {
    CpvAccess(Node_instance) = this;
    eventEndOfTimeStep = traceRegisterUserEvent("EndOfTimeStep");
  } else {
    NAMD_bug("Node::Node() - another instance of Node exists!");
  }

  CpvAccess(BOCclass_group) = msg->group;
  delete msg;

  CpvAccess(BOCclass_group).node = thisgroup;

  startupPhase = 0;
  numNodeStartup = CkNumPes();

  molecule = NULL;
  parameters = NULL;
  simParameters = NULL;
  configList = NULL;
  pdb = NULL;
  state = NULL;
  output = NULL;
  imd = new IMDOutput;

  Compute::setNode(this);

  DebugM(4,"Creating PatchMap, AtomMap, ComputeMap\n");
  patchMap = PatchMap::Instance();
  atomMap = AtomMap::Instance();
  computeMap = ComputeMap::Instance();

  DebugM(4,"Binding to BOC's\n");
  CProxy_PatchMgr pm(CpvAccess(BOCclass_group).patchMgr);
  patchMgr = pm.ckLocalBranch();
  CProxy_ProxyMgr prm(CpvAccess(BOCclass_group).proxyMgr);
  proxyMgr = prm.ckLocalBranch();
  CProxy_WorkDistrib wd(CpvAccess(BOCclass_group).workDistrib);
  workDistrib = wd.ckLocalBranch();
  CProxy_ComputeMgr cm(CpvAccess(BOCclass_group).computeMgr);
  computeMgr = cm.ckLocalBranch();
  CProxy_LdbCoordinator lc(CpvAccess(BOCclass_group).ldbCoordinator);
  ldbCoordinator = lc.ckLocalBranch();

// #if !defined(NOHOSTNAME)
//   // Where are we?
//   char host[1024];
//   gethostname(host, 1024);
//   iout << iINFO << iPE << " Starting out on host: " << host << "\n" << endi;
// #endif
}

//----------------------------------------------------------------------
// ~Node(void) needs to clean up everything.

Node::~Node(void)
{
  delete output;
  delete computeMap;
  delete atomMap;
  delete patchMap;
  delete CpvAccess(comm);
}

//----------------------------------------------------------------------
// Startup Sequence

void Node::messageStartUp() {
  (CProxy_Node(CpvAccess(BOCclass_group).node)).startup();
}

void Node::startUp(CkQdMsg *qmsg) {
  delete qmsg;
  (CProxy_Node(CpvAccess(BOCclass_group).node)).startup();
}


void Node::startup() {
  int gotoRun = false;

  if (!CkMyPe()) {
     iout << iINFO << "Entering startup phase " << startupPhase << " with " <<
	(memusage()/1024) << " kB of memory in use.\n" << endi;
  }
  
  switch (startupPhase) {

  case 0:
    namdOneCommInit(); // Namd1.X style
  break;

  case 1:
    // send & receive molecule, simparameters... (Namd1.X style)
    if (CkMyPe()) {
      namdOneRecv();
    } else {
      namdOneSend();
    }
  break;

  case 2:
    // take care of inital thread setting
    threadInit();

    // create blank AtomMap
    AtomMap::Object()->allocateMap(molecule->numAtoms);
  break;

  case 3:
    if (!CkMyPe()) {
      output = new Output; // create output object just on PE(0)
      workDistrib->patchMapInit(); // create space division
      workDistrib->createHomePatches(); // load atoms into HomePatch(es)
      workDistrib->assignNodeToPatch();
      workDistrib->mapComputes();
      ComputeMap::Object()->printComputeMap();
      workDistrib->sendMaps();
    }
  break;

  case 4:
    if (!CkMyPe()) {
      workDistrib->distributeHomePatches();
    }
  break;

  case 5: 
    Sync::Object()->openSync();  // decide if to open local Sync 
    proxyMgr->createProxies();  // need Home patches before this
  break;

  case 6:
    if (!CkMyPe()) {
      ComputeMap::Object()->printComputeMap();
    }
    if (proxySendSpanning || proxyRecvSpanning )
      proxyMgr->buildProxySpanningTree();
    DebugM(4,"Creating Computes\n");
    computeMgr->createComputes(ComputeMap::Object());
    DebugM(4,"Building Sequencers\n");
    buildSequencers();
    DebugM(4,"Initializing LDB\n");
    LdbCoordinator::Object()->initialize(patchMap,computeMap);
  break;

  case 7:
    gotoRun = true;
  break;

  default:
    NAMD_bug("Startup Phase has a bug - check case statement");
  break;

  }

  startupPhase++;
  if (!CkMyPe()) {
    if (!gotoRun) {
      CkStartQD(CProxy_Node::ckIdx_startUp((CkQdMsg*)0),&thishandle);
    } else {
      Node::messageRun();
    }
  }
}

void Node::namdOneCommInit()
{
  if (CpvAccess(comm) == NULL) {
    CpvAccess(comm) = new Communicate();
#ifdef DPMTA
    pvmc_init();
#endif
  }
}

// Namd 1.X style Send/Recv of simulation information

void Node::namdOneRecv() {
  MIStream *conv_msg;

  // Receive molecule and simulation parameter information
  simParameters = new SimParameters;
  //****** BEGIN CHARMM/XPLOR type changes
  parameters = new Parameters(simParameters);
  //****** END CHARMM/XPLOR type changes
  molecule = new Molecule(simParameters,parameters);

  DebugM(4, "Getting SimParameters\n");
  conv_msg = CpvAccess(comm)->newInputStream(0, SIMPARAMSTAG);
  simParameters->receive_SimParameters(conv_msg);

  DebugM(4, "Getting Parameters\n");
  conv_msg = CpvAccess(comm)->newInputStream(0, STATICPARAMSTAG);
  parameters->receive_Parameters(conv_msg);

  DebugM(4, "Getting Molecule\n");
  conv_msg = CpvAccess(comm)->newInputStream(0, MOLECULETAG);
  molecule->receive_Molecule(conv_msg);

  DebugM(4, "Done Receiving\n");
}

void Node::namdOneSend() {
  // I'm Pe(0) so I send what I know
  DebugM(4, "Sending SimParameters\n");
  simParameters->send_SimParameters(CpvAccess(comm));
  DebugM(4, "Sending Parameters\n");
  parameters->send_Parameters(CpvAccess(comm));
  DebugM(4, "Sending Molecule\n");
  molecule->send_Molecule(CpvAccess(comm));
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
  (CProxy_Node(CpvAccess(BOCclass_group).node)).run();
}


//-----------------------------------------------------------------------
// run(void) runs the specified simulation for the specified number of
// steps, overriding the contents of the configuration file
//-----------------------------------------------------------------------
void Node::run()
{
  // iout << iINFO << iPE << " Running\n" << endi;
  numHomePatchesRunning = patchMap->numHomePatches();
  if (CkMyPe() == 0) {
    numHomePatchesRunning++; //Take into account controller on node 0
  } 

  //Check if the number of patches is less than the number of nodes
  if (numNodes() > patchMap->numPatches())
	numNodesRunning=patchMap->numPatches();
  else 
  	numNodesRunning = numNodes();

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
     iout << iINFO << "Finished startup with " <<
	(memusage()/1024) << " kB of memory in use.\n" << endi;
  }
  
}


//-----------------------------------------------------------------------
// Node scriptBarrier() - twiddle parameters with simulation halted
//-----------------------------------------------------------------------

void Node::enableScriptBarrier() {
  CkStartQD(CProxy_Node::ckIdx_scriptBarrier((CkQdMsg*)0),&thishandle);
}

void Node::scriptBarrier(CkQdMsg *qmsg) {
  delete qmsg;
  //script->awaken();
}

void Node::scriptParam(ScriptParamMsg *msg) {
  simParameters->scriptSet(msg->param,msg->value);
  delete msg;
}


void Node::sendEnableExitScheduler(void) {
  //CmiPrintf("sendEnableExitScheduler\n");
  CkQdMsg *msg = new CkQdMsg;
  CProxy_Node(thisgroup).recvEnableExitScheduler(msg,0);
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
    CkStartQD(CProxy_Node::ckIdx_exitScheduler((CkQdMsg*)0),&thishandle);
  }
}

void Node::exitScheduler(CkQdMsg *msg) {
  //CmiPrintf("exitScheduler %d\n",CkMyPe());
  CsdExitScheduler();
  delete msg;
}

void Node::sendEnableEarlyExit(void) {
  CkQdMsg *msg = new CkQdMsg;
  CProxy_Node(thisgroup).recvEnableEarlyExit(msg,0);
}

void Node::recvEnableEarlyExit(CkQdMsg *msg) {
  delete msg;
  enableEarlyExit();
}

void Node::enableEarlyExit(void) {
  if ( CkMyPe() ) {
    sendEnableEarlyExit();
  } else {
    CkStartQD(CProxy_Node::ckIdx_earlyExit((CkQdMsg*)0),&thishandle);
  }
}

void Node::earlyExit(CkQdMsg *msg) {
  iout << iERROR << "Exiting prematurely.\n" << endi;
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

//======================================================================
// Private functions


#include "Node.def.h"

