/***************************************************************************/
/*          (C) Copyright 1996,1997 The Board of Trustees of the           */
/*                          University of Illinois                         */
/*                           All Rights Reserved                           */
/***************************************************************************/
/***************************************************************************
 * DESCRIPTION: Toplevel routines for initializing a Node for a simulation
 *              one Node per Pe (processor element).
 *
 ***************************************************************************/

static char ident[] = "@(#)$Header: /home/cvs/namd/cvsroot/namd2/src/Node.C,v 1.1015 1997/04/06 22:45:06 ari Exp $";

#include <unistd.h>
#include "ckdefs.h"
#include "chare.h"
#include "c++interface.h"
#include "Node.top.h"
#include "Node.h"
#include "Namd.h"

#define MIN_DEBUG_LEVEL 3
//#define DEBUGM
#include "Debug.h"

#include <stdio.h>
#include "converse.h"

#include "unistd.h"

#include "main.top.h"
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
#include "Parameters.h"
#include "SimParameters.h"
#include "CommunicateConverse.h"
#include "Inform.h"
#include "LdbCoordinator.h"

extern Communicate *comm;

BOCgroup BOCclass::group;

//======================================================================
// Public Functions

//----------------------------------------------------------------------
// Singleton implementation
Node *Node::_instance = 0;

//----------------------------------------------------------------------
// BOC constructor
Node::Node(GroupInitMsg *msg)
{
  if (_instance == 0) {
    _instance = this;
  } else {
    iout << iERRORF << "Node::Node() - another instance of Node exists!\n"
      << endi;
    Namd::die();
  }

  group = msg->group;
  delete msg;

  group.node = thisgroup;

  startupPhase = 0;
  numNodeStartup = CNumPes();

  molecule = NULL;
  parameters = NULL;
  simParameters = NULL;
  configList = NULL;
  pdb = NULL;
  state = NULL;
  output = NULL;

  Compute::setNode(this);

  patchMap = PatchMap::Instance();
  atomMap = AtomMap::Instance();
  computeMap = ComputeMap::Instance();

  patchMgr = CLocalBranch(PatchMgr,group.patchMgr);
  proxyMgr = CLocalBranch(ProxyMgr,group.proxyMgr);
  workDistrib = CLocalBranch(WorkDistrib,group.workDistrib);
  computeMgr = CLocalBranch(ComputeMgr,group.computeMgr);
  ldbCoordinator = CLocalBranch(LdbCoordinator,group.ldbCoordinator);

  // Where are we?
  char host[1024];
  gethostname(host, 1024);
  iout << iINFO << iPE << " Starting out on host: " << host << "\n" << endi;
}

//----------------------------------------------------------------------
// ~Node(void) needs to clean up everything.

Node::~Node(void)
{
  delete output;
  delete computeMap;
  delete atomMap;
  delete patchMap;
  delete comm;
}

//----------------------------------------------------------------------
// Startup Sequence

void Node::messageStartUp() {
  InitMsg *msg = new (MsgIndex(InitMsg)) InitMsg;
  CBroadcastMsgBranch(Node, startup, msg, group.node);
}

void Node::startUp(QuiescenceMessage *qmsg) {
  delete qmsg;
  InitMsg *msg = new (MsgIndex(InitMsg)) InitMsg;
  CBroadcastMsgBranch(Node, startup, msg, group.node);
}


void Node::startup(InitMsg *msg) {
  delete msg;

  int gotoRun = false;
  
  switch (startupPhase) {

  case 0:
    if (!CMyPe()) {
       iout << iINFO << " Starting Phase " << startupPhase << "\n" << endi;
    }
    namdOneCommInit(); // Namd1.X style
  break;

  case 1:
    if (!CMyPe()) {
       iout << iINFO << " Starting Phase " << startupPhase << "\n" << endi;
    }

    // send & receive molecule, simparameters... (Namd1.X style)
    if (CMyPe()) {
      namdOneRecv();
    } else {
      namdOneSend();
    }
  break;

  case 2:
    if (!CMyPe()) {
       iout << iINFO << " Starting Phase " << startupPhase << "\n" << endi;
    }

    // take care of inital thread setting
    threadInit();

    // create blank AtomMap
    AtomMap::Object()->allocateMap(molecule->numAtoms);
  break;

  case 3:
    if (!CMyPe()) {
       iout << iINFO << " Starting Phase " << startupPhase << "\n" << endi;
    }
    if (!CMyPe()) {
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
    if (!CMyPe()) {
       iout << iINFO << " Starting Phase " << startupPhase << "\n" << endi;

      workDistrib->distributeHomePatches();
    }
  break;

  case 5: 
    if (!CMyPe()) {
       iout << iINFO << " Starting Phase " << startupPhase << "\n" << endi;
    }
    proxyMgr->createProxies();  // need Home patches before this
  break;

  case 6:
    if (!CMyPe()) {
       iout << iINFO << " Starting Phase " << startupPhase << "\n" << endi;
      ComputeMap::Object()->printComputeMap();
    }
    DebugM(4,"Creating Computes\n");
    computeMgr->createComputes(ComputeMap::Object());
    DebugM(4,"Building Sequencers\n");
    buildSequencers();
    DebugM(4,"Initializing LDB\n");
    LdbCoordinator::Object()->initialize(patchMap,computeMap);
  break;

  case 7:
    if (!CMyPe()) {
       iout << iINFO << " Starting Phase " << startupPhase << "\n" << endi;
    }
    gotoRun = true;
  break;

  default:
    iout << iERRORF << iPE << "Startup Phase has a bug - check case statement\n" << endi;
    Namd::die();
  break;

  }

  startupPhase++;
  if (!CMyPe()) {
    if (!gotoRun) {
      CStartQuiescence(GetEntryPtr(Node,startUp), thishandle);
    } else {
      Node::messageRun();
      // CStartQuiescence(GetEntryPtr(Node,quiescence), thishandle);
    }
  }
}

void Node::namdOneCommInit()
{
  if (comm == NULL) {
    comm = new CommunicateConverse(0,0,1);
  }
}

// Namd 1.X style Send/Recv of simulation information

void Node::namdOneRecv() {
  Message *conv_msg=NULL;
  int tag;
  int zero=0;

  // Receive molecule and simulation parameter information
  simParameters = new SimParameters;
  parameters = new Parameters;
  molecule = new Molecule(simParameters);

  DebugM(4, "Getting SimParameters\n");
  do{
    tag=SIMPARAMSTAG;
    conv_msg = comm->receive(zero,tag);
  } while (conv_msg == NULL);
  simParameters->receive_SimParameters(conv_msg);

  DebugM(4, "Getting Parameters\n");
  do{
    tag=STATICPARAMSTAG;
    conv_msg = comm->receive(zero,tag);
  } while (conv_msg == NULL);
  parameters->receive_Parameters(conv_msg);

  DebugM(4, "Getting Molecule\n");
  do{
    tag=MOLECULETAG;
    conv_msg = comm->receive(zero,tag);
  } while (conv_msg == NULL);
  molecule->receive_Molecule(conv_msg);

  DebugM(4, "Done Receiving\n");
}

void Node::namdOneSend() {
  // I'm Pe(0) so I send what I know
  DebugM(4, "Sending SimParameters\n");
  simParameters->send_SimParameters(comm);
  DebugM(4, "Sending Parameters\n");
  parameters->send_Parameters(comm);
  DebugM(4, "Sending Molecule\n");
  molecule->send_Molecule(comm);
  DebugM(4, "Done Sending\n");
}

// Initial thread setup

void Node::threadInit() {
  // Thread initialization
  if (CthImplemented()) {
    CthSetStrategyDefault(CthSelf());
  } else {
    iout << iERRORF 
      << "Node::startup() Oh no, tiny elvis, threads not implemented\n"
      << endi;
    Namd::die();
  }
}

//
void Node::buildSequencers() {
  HomePatchList *hpl = PatchMap::Object()->homePatchList();
  ResizeArrayIter<HomePatchElem> ai(*hpl);

  // Controller object is only on Pe(0)
  if ( ! CMyPe() ) {
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
  RunMsg *msg = new (MsgIndex(RunMsg)) RunMsg;
  CBroadcastMsgBranch(Node, run, msg, group.node);
}


//-----------------------------------------------------------------------
// run(void) runs the specified simulation for the specified number of
// steps, overriding the contents of the configuration file
//-----------------------------------------------------------------------
void Node::run(RunMsg *msg)
{
  delete msg;

  iout << iINFO << iPE << " Running\n" << endi;
  numHomePatchesRunning = patchMap->numHomePatches();
  if (CMyPe() == 0) {
    numHomePatchesRunning++; //Take into account controller on node 0
  } 
  numNodesRunning = numNodes();

  Namd::startTimer();  // We count timings from this point on

  // Start Controller (aka scalar Sequencer) on Pe(0)
  if ( ! CMyPe() ) {
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
}


//-----------------------------------------------------------------------
// Node homeDone() - broadcast to all nodes
//-----------------------------------------------------------------------
void Node::messageHomeDone() {
  DoneMsg *msg = new (MsgIndex(DoneMsg)) DoneMsg;
  CSendMsgBranch(Node, homeDone, msg, group.node, CMyPe());
}


//-----------------------------------------------------------------------
// per node countdown of HomePatch completions.  Signals nodeDone()
// when completed.
//-----------------------------------------------------------------------
void Node::homeDone(DoneMsg *msg) {
  delete msg;

  if (!--numHomePatchesRunning) {
     DoneMsg *msg = new (MsgIndex(DoneMsg)) DoneMsg;
     CSendMsgBranch(Node, nodeDone, msg, group.node, 0);
  }
}


//-----------------------------------------------------------------------
// Countdown of Node completions - 
//-----------------------------------------------------------------------
void Node::nodeDone(DoneMsg *msg) {
  delete msg;

  if (!--numNodesRunning && !CMyPe()) {
     Namd::namdDone();
  }
}


//-----------------------------------------------------------------------
// Deal with quiescence - this terminates the program (for now)
//-----------------------------------------------------------------------
void Node::quiescence(QuiescenceMessage * msg)
{
  delete msg;

  iout << iINFO << iPE << "Quiescence detected, exiting Charm.\n" << endi;
  CharmExit();
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


#include "Node.bot.h"

/***************************************************************************
 * RCS INFORMATION:
 *
 *	$RCSfile: Node.C,v $
 *	$Author: ari $	$Locker:  $		$State: Exp $
 *	$Revision: 1.1015 $	$Date: 1997/04/06 22:45:06 $
 *
 ***************************************************************************
 * REVISION HISTORY:
 *
 * $Log: Node.C,v $
 * Revision 1.1015  1997/04/06 22:45:06  ari
 * Add priorities to messages.  Mods to help proxies without computes.
 * Added quick enhancement to end of list insertion of ResizeArray(s)
 *
 * Revision 1.1014  1997/04/04 23:34:22  milind
 * Got NAMD2 to run on Origin2000.
 * Included definitions of class static variables in C files.
 * Fixed alignment bugs by using memcpy instead of assignment in
 * pack and unpack.
 *
 * Revision 1.1013  1997/04/04 22:28:54  ari
 * Removed quiescence detection after the initial startup phase.
 *
 * Revision 1.1012  1997/03/27 20:25:49  brunner
 * Changes for LdbCoordinator, the load balance control BOC
 *
 * Revision 1.1011  1997/03/21 16:36:00  nealk
 * Modified sorting case when both atoms are group members.
 * Was incorrectly sorting them by their atomID.  Now sorting by GPID.
 * Also turned off debugging in Node.C.
 *
 * Revision 1.1010  1997/03/20 23:53:45  ari
 * Some changes for comments. Copyright date additions.
 * Hooks for base level update of Compute objects from ComputeMap
 * by ComputeMgr.  Useful for new compute migration functionality.
 *
 * Revision 1.1009  1997/03/19 11:54:35  ari
 * Add Broadcast mechanism.
 * Fixed RCS Log entries on files that did not have Log entries.
 * Added some register variables to Molecule and ComputeNonbondedExcl.C
 *
 * Revision 1.1008  1997/03/14 21:40:11  ari
 * Reorganized startup to make possible inital load
 * balancing by changing methods in WorkDistrib.
 * Also made startup more transparent and easier
 * to modify.
 *
 * Revision 1.1007  1997/03/10 17:40:14  ari
 * UniqueSet changes - some more commenting and cleanup
 *
 * Revision 1.1006  1997/03/06 22:06:06  ari
 * Removed Compute.ci
 * Comments added - more code cleaning
 *
 * Revision 1.1005  1997/03/04 22:37:15  ari
 * Clean up of code.  Debug statements removal, dead code removal.
 * Minor fixes, output fixes.
 * Commented some code from the top->down.  Mainly reworked Namd, Node, main.
 *
 * Revision 1.1004  1997/02/26 16:53:12  ari
 * Cleaning and debuging for memory leaks.
 * Adding comments.
 * Removed some dead code due to use of Quiescense detection.
 *
 * Revision 1.1003  1997/02/13 16:17:16  ari
 * Intermediate debuging commit - working to fix deep bug in migration?
 *
 * Revision 1.1002  1997/02/13 04:43:10  jim
 * Fixed initial hanging (bug in PatchMap, but it still shouldn't have
 * happened) and saved migration messages in the buffer from being
 * deleted, but migration still dies (even on one node).
 *
 * Revision 1.1001  1997/02/11 22:56:13  jim
 * Added dcd file writing.
 *
 * Revision 1.1000  1997/02/06 15:58:50  ari
 * Resetting CVS to merge branches back into the main trunk.
 * We will stick to main trunk development as suggested by CVS manual.
 * We will set up tags to track fixed points of development/release
 * as suggested by CVS manual - all praise the CVS manual.
 *
 * Revision 1.779  1997/02/06 15:53:18  ari
 * Updating Revision Line, getting rid of branches
 *
 * Revision 1.778.2.1  1997/02/05 22:18:15  ari
 * Added migration code - Currently the framework is
 * there with compiling code.  This version does
 * crash shortly after migration is complete.
 * Migration appears to complete, but Patches do
 * not appear to be left in a correct state.
 *
 * Revision 1.778  1997/01/28 00:30:59  ari
 * internal release uplevel to 1.778
 *
 * Revision 1.777.2.1  1997/01/19 21:45:01  jim
 * Removed some debug messages.
 *
 * Revision 1.777  1997/01/17 19:36:34  ari
 * Internal CVS leveling release.  Start development code work
 * at 1.777.1.1.
 *
 * Revision 1.35  1997/01/10 21:33:22  milind
 * Added some debugging statements
 *
 * Revision 1.34  1997/01/10 17:41:25  jim
 * *** empty log message ***
 *
 * Revision 1.33  1997/01/09 20:48:10  jim
 * added Controller code
 *
 * Revision 1.32  1996/12/30 20:37:38  nealk
 * Mode debug code.
 *
 * Revision 1.31  1996/12/27 22:22:33  nealk
 * Made extern comm global to file (was already) and readded second comm=new.
 *
 * Revision 1.30  1996/12/26 22:26:54  nealk
 * Corrected one of the communications bugs for +p2 -- comm wasn't initialized
 * at the right time.
 *
 * Revision 1.29  1996/12/23 17:37:26  nealk
 * Mode debug code...
 *
 * Revision 1.28  1996/12/20 22:53:27  jim
 * fixing parallel bugs, going home for Christmass
 *
 * Revision 1.27  1996/12/19 02:26:30  jim
 * Node::startup2 is now triggered by quiescence
 *
 * Revision 1.26  1996/12/16 19:06:30  nealk
 * Debugging a core dump when it registers force functions.
 *
 * Revision 1.25  1996/12/13 08:52:37  jim
 * staged startup implemented
 *
 * Revision 1.24  1996/12/12 20:14:50  milind
 * *** empty log message ***
 *
 * Revision 1.23  1996/12/12 08:58:12  jim
 * startup corrections (I hope).
 *
 * Revision 1.22  1996/12/11 22:24:13  jim
 * made node startup and workDistrib calls work in parallel
 *
 * Revision 1.21  1996/12/11 00:04:23  milind
 * *** empty log message ***
 *
 * Revision 1.20  1996/12/10 00:13:12  ari
 * *** empty log message ***
 *
 * Revision 1.19  1996/12/06 19:54:12  ari
 * *** empty log message ***
 *
 * Revision 1.18  1996/12/02 19:39:26  nealk
 * Removed DEBUGM macro
 *
 * Revision 1.17  1996/12/01 22:46:11  jim
 * switched to use simParams for number of cycles
 *
 * Revision 1.16  1996/11/30 00:44:24  jim
 * added sequencer use, ComputeMgr use, and quiescence detection
 *
 * Revision 1.15  1996/11/26 16:33:35  nealk
 * Added ComputeDihedrals and ComputeImpropers.
 *
 * Revision 1.14  1996/11/23 23:01:00  jim
 * added ComputeNonbondedSelf test code
 *
 * Revision 1.13  1996/11/22 01:44:53  jim
 * added calls to service AtomMap
 *
 * Revision 1.12  1996/11/22 01:02:18  ari
 * *** empty log message ***
 *
 * Revision 1.11  1996/11/22 00:18:51  ari
 * *** empty log message ***
 *
 * Revision 1.10  1996/11/05 16:59:58  ari
 * *** empty log message ***
 *
 * Revision 1.9  1996/10/29 23:35:27  ari
 * *** empty log message ***
 *
 * Revision 1.8  1996/10/29 17:18:34  brunner
 * Changed workDistrib calls
 *
 * Revision 1.7  1996/10/16 08:22:39  ari
 * *** empty log message ***
 *
 * Revision 1.6  1996/10/04 22:24:15  brunner
 * Took out call to createPatches.  Maybe this will have to go back in
 *
 * Revision 1.5  1996/08/23 22:03:52  brunner
 * Made WorkdDistrib, PatchMgr public members
 *
 * Revision 1.4  1996/08/21 23:58:25  brunner
 * *** empty log message ***
 *
 * Revision 1.3  1996/08/19 22:05:31  ari
 * *** empty log message ***
 *
 * Revision 1.2  1996/08/16 21:42:29  brunner
 * *** empty log message ***
 *
 * Revision 1.1  1996/08/16 21:19:34  ari
 * Initial revision
 *
 ***************************************************************************/
