/***************************************************************************/
/*                                                                         */
/*              (C) Copyright 1996 The Board of Trustees of the            */
/*                          University of Illinois                         */
/*                           All Rights Reserved                           */
/*                                                                         */
/***************************************************************************/

/***************************************************************************
 * DESCRIPTION:
 *
 ***************************************************************************/

static char ident[] = "@(#)$Header: /home/cvs/namd/cvsroot/namd2/src/Node.C,v 1.1002 1997/02/13 04:43:10 jim Exp $";


#include "ckdefs.h"
#include "chare.h"
#include "c++interface.h"
#include "Node.top.h"
#include "Node.h"

#include <stdio.h>
#include "converse.h"

#include "unistd.h"

#include "main.top.h"
#include "main.h"
#include "WorkDistrib.h"
#include "PatchMgr.h"
#include "Patch.h"
#include "Compute.h"
#include "ComputeAngles.h"
#include "ComputeDihedrals.h"
#include "ComputeImpropers.h"
#include "ComputeNonbondedSelf.h"
#include "ComputeMap.h"
#include "ComputeMgr.h"
#include "Molecule.h"
#include "HomePatchList.h"
#include "AtomMap.h"
#include "Sequencer.h"
#include "Controller.h"
#include "NamdState.h"
#include "Output.h"
//#include "ProxyMgr.h"
//#include "MessageComm.h"
//#include "PatchMap.h"
#include "Parameters.h"
#include "SimParameters.h"
#include "CommunicateConverse.h"

extern Communicate *comm;

#define MIN_DEBUG_LEVEL 3
#define DEBUGM
#include "Debug.h"

//======================================================================
// Public Functions


Node *Node::_instance = 0;

Node::Node(GroupInitMsg *msg)
{
  DebugM(1, "Node::Node() - starting\n");
  group = msg->group;
  group.node = thisgroup;
  delete msg;

  /* in case any node other than 0 reaches this place first. */
  if (comm == NULL)
  {
    comm = new CommunicateConverse(0,0);
    DebugM(3,"comm being initialized.\n");
  }
  else DebugM(3,"comm already initialized.\n");

  molecule = NULL;
  parameters = NULL;
  simParameters = NULL;
  configList = NULL;
  pdb = NULL;
  state = NULL;
  output = NULL;
  Compute::setNode(this);

  if (_instance == 0) {
    _instance = this;
  } else {
    DebugM(1, "Node::Node() - another instance of Node exists!\n");
  }

  DebugM(1, "Node::Node() - constructing\n");

  patchMap = PatchMap::Instance();
  atomMap = AtomMap::Instance();
  computeMap = ComputeMap::Instance();

  DebugM(1,"Node::Node() - I have the following groups\n");
  DebugM(1,"Node::Node() - WorkDistrib " << group.workDistrib << "\n");
  DebugM(1,"Node::Node() - PatchMgr " << group.patchMgr << "\n");

  // Put num<name of BOC>Startup here
  numNodeStartup = CNumPes();
#if 0
  Message *conv_msg=NULL;
  conv_msg = new Message;
  comm->send_method(Communicate::NOW);
  conv_msg->put(1);
  if (CMyPe()) {
    DebugM(1,"Sending checkin to node 0\n");
    int rc = comm->send(conv_msg,0,DISTRIBTAG);
    DebugM(1,"Sent checkin to node 0, rc=" << rc << "\n");
  }
#endif
}

//----------------------------------------------------------------------
// ~Node(void) needs to clean up everything.

Node::~Node(void)
{
  delete output;
}

//----------------------------------------------------------------------
// Node startup
//
// startup(void) receives all the startup data from the Namd object,
// and does everything necessary to run the simulation
//    Ask work_distrib for initial distribution
//    This must be split-phase, so results can be received from node 0.
// Receipt of the maps send by sendMaps() triggers the patch and compute
// object creation

void Node::messageStartup() {
  InitMsg *msg = new (MsgIndex(InitMsg)) InitMsg;
  DebugM(1,"calling CBroadcastMsgBranch()\n");
  CBroadcastMsgBranch(Node, startup, msg, group.node);
  DebugM(1,"returned from CBroadcastMsgBranch()\n");
}

void Node::startup(InitMsg *msg)
{
  char **argvdummy;
  Message *conv_msg=NULL;
  
  delete msg;

  /* in case node 0 reaches here before other nodes... */
  if (comm == NULL)
  {
    comm = new CommunicateConverse(0,0);
    DebugM(3,"comm being initialized.\n");
  }
  else DebugM(3,"comm already initialized.\n");

  if ( CMyPe() ) {
     simParameters = new SimParameters;
     parameters = new Parameters;
     molecule = new Molecule(simParameters);
     do{
        conv_msg = comm->receive(0,SIMPARAMSTAG);
     } while (conv_msg == NULL);
     simParameters->receive_SimParameters(conv_msg);
     do{
        conv_msg = comm->receive(0,STATICPARAMSTAG);
     } while (conv_msg == NULL);
     parameters->receive_Parameters(conv_msg);
     do{
        conv_msg = comm->receive(0,MOLECULETAG);
     } while (conv_msg == NULL);
     molecule->receive_Molecule(conv_msg);
  } else {
     simParameters->send_SimParameters(comm);
     parameters->send_Parameters(comm);
     molecule->send_Molecule(comm);
  }

  // Thread initialization
  if (CthImplemented()) {
    CthSetStrategyDefault(CthSelf());
  } else {
    DebugM(10, "Node::startup() - Oh no, tiny elvis, threads not implemented\n");
    CharmExit();
  }

  AtomMap::Object()->allocateMap(molecule->numAtoms);

  DebugM(1,iPE << " Calling CLocalBranch(WorkDistrib)\n");
  workDistrib = CLocalBranch(WorkDistrib,group.workDistrib);

  if ( CMyPe() ) {
     conv_msg = new Message;
     conv_msg->put(1);
     comm->send_method(Communicate::NOW);
     DebugM(1,"Sending checkin to node 0 from node " << iPE << "\n");
     comm->send(conv_msg,0,DISTRIBTAG);
     DebugM(1,"Sent checkin to node 0 from node " << iPE << "\n");
  } else {
     for ( int i = 1; i < numNodes(); ++i )
     {
       DebugM(1,"Looking for checkin from node " << i << "\n");
       do
       {
         conv_msg = comm->receive(i,DISTRIBTAG);
       } while (conv_msg == NULL);
       delete conv_msg;
       DebugM(1,"Received checkin from node " << i << "\n");
     }
     DebugM(1,"Received all checkins, messaging startup1()\n");
     InitMsg *msg = new (MsgIndex(InitMsg)) InitMsg;
     CSendMsgBranch(Node, startup1, msg, group.node, 0);
  }
}

void Node::startup1(InitMsg *msg)
{
  DebugM(1,"In startup1() on node " << CMyPe() << endl);
  delete msg;

  DebugM(1, "workDistrib->buildMaps() Pe=" << CMyPe() << "\n");
  workDistrib->buildMaps();
  DebugM(1, "workDistrib->sendMaps() Pe=" << CMyPe() << "\n");
  workDistrib->sendMaps();

  DebugM(1,"CStartQuiescence in startup1\n");
  CStartQuiescence(GetEntryPtr(Node,messageStartup2), thishandle);
  DebugM(1,"End of startup1()\n");
}

void Node::messageStartup2(QuiescenceMessage * qm) {
  DebugM(1,"In messageStartup2() on node " << CMyPe() << endl);
  delete qm;
  InitMsg *msg = new (MsgIndex(InitMsg)) InitMsg;
  CBroadcastMsgBranch(Node, startup2, msg, group.node);
}

void Node::startup2(InitMsg *msg)
{
  DebugM(1,"In startup2() on node " << CMyPe() << endl);
  delete msg;
  
  DebugM(1, "workDistrib->createPatches() Pe=" << CMyPe() << "\n");
  if ( ! CMyPe() )
  {
    ComputeMap::Object()->printComputeMap();
    workDistrib->createPatches();
    output = new Output;
  }

  if ( ! CMyPe() )
  {
	DebugM(1,"CStartQuiescence in startup2\n");
	CStartQuiescence(GetEntryPtr(Node,messageStartup3), thishandle);
  }
  DebugM(1,"End of startup2()\n");
}

void Node::messageStartup3(QuiescenceMessage * qm) {
  DebugM(1,"In messageStartup3() on node " << CMyPe() << endl);
  delete qm;
  InitMsg *msg = new (MsgIndex(InitMsg)) InitMsg;
  CBroadcastMsgBranch(Node, startup3, msg, group.node);
}

void Node::startup3(InitMsg *msg)
{
  DebugM(1,"In startup3() on node " << CMyPe() << endl);
  delete msg;

  // create proxies
  proxyMgr = CLocalBranch(ProxyMgr,group.proxyMgr);
  proxyMgr->createProxies();

  if ( ! CMyPe() )
  {
	DebugM(1,"CStartQuiescence in startup3\n");
	CStartQuiescence(GetEntryPtr(Node,messageStartup4), thishandle);
  }
  DebugM(1,"End of startup3()\n");
}

void Node::messageStartup4(QuiescenceMessage * qm) {
  DebugM(1,"In messageStartup4() on node " << CMyPe() << endl);
  delete qm;
  InitMsg *msg = new (MsgIndex(InitMsg)) InitMsg;
  CBroadcastMsgBranch(Node, startup4, msg, group.node);
}

void Node::startup4(InitMsg *msg)
{
  DebugM(1,"In startup4() on node " << CMyPe() << endl);
  delete msg;

  // create computes
  computeMgr = CLocalBranch(ComputeMgr,group.computeMgr);
  DebugM(3, "Trying to create computes.\n");
  computeMgr->createComputes(ComputeMap::Object());
  DebugM(3, "Created computes.  Making patch managers\n");
  patchMgr = CLocalBranch(PatchMgr,group.patchMgr);
  DebugM(3, "Created patch managers\n");

  HomePatchList *hpl = PatchMap::Object()->homePatchList();
  ResizeArrayIter<HomePatchElem> ai(*hpl);

  if ( ! CMyPe() )
  {
    DebugM(4, "Node::startup4() - creating controller.\n");
    Controller *controller = new Controller(state);
    state->useController(controller);
  }

  DebugM(4, "Node::startup4() - iterating over home patches!\n");
  for (ai=ai.begin(); ai != ai.end(); ai++) {
    HomePatch *p = (*ai).p;
    DebugM(1, "Node::run() - signaling patch "<< p->getPatchID() << endl);
    Sequencer *sequencer = new Sequencer(p);
    p->useSequencer(sequencer);
  }

  messageStartupDone();   // collect on master node
  DebugM(1,"End of startup4()\n");
}


//-----------------------------------------------------------------------
// Node Startup completion messaging and barrier on node 0

void Node::messageStartupDone() {
  DoneMsg *msg = new (MsgIndex(DoneMsg)) DoneMsg;
  CSendMsgBranch(Node, startupDone, msg, group.node, 0);
}

void Node::startupDone(DoneMsg *msg) {
  delete msg;

  if (CMyPe() == 0) {
    DebugM(1, "Node::startupDone() - got one\n");
    if (!--numNodeStartup) {
      DebugM(1, "Node::startupDone() - triggered run() \n");
      Node::messageRun();
      DebugM(1,"CStartQuiescence in startupDone\n");
      CStartQuiescence(GetEntryPtr(Node,quiescence), thishandle);
    }
  } else {
    DebugM(1, "Node::startupDone() - message sent to wrong Pe!\n");
  }
}

//-----------------------------------------------------------------------
// Node run() - broadcast to all nodes

void Node::messageRun() {
  RunMsg *msg = new (MsgIndex(RunMsg)) RunMsg;
  CBroadcastMsgBranch(Node, run, msg, group.node);
}


// run(void) runs the specified simulation for the specified number of
// steps, overriding the contents of the configuration file
void Node::run(RunMsg *msg)
{
  delete msg;

  // This is testbed code!
  DebugM(4, "Node::run() - invoked\n");
  DebugM(1, "Node::run() - message address was " << msg << "\n");
  static int foo = 0;
  if ( ! foo ) foo = 1;
  else return;

  HomePatchList *hpl = PatchMap::Object()->homePatchList();
  ResizeArrayIter<HomePatchElem> ai(*hpl);

  if ( ! CMyPe() )
  {
    state->runController();
  }

  for (ai=ai.begin(); ai != ai.end(); ai++) {
    HomePatch *p = (*ai).p;
    p->runSequencer();
  }
}


// Deal with quiescence - this terminates the program (for now)
void Node::quiescence(QuiescenceMessage * msg)
{
  delete msg;

  DebugM(4, "Quiescence detected, exiting Charm.\n");
  CharmExit();
}


//------------------------------------------------------------------------
// Some odd utilities

void Node::saveMolDataPointers(Molecule *molecule,
			       Parameters *parameters,
			       SimParameters *simParameters,
			       ConfigList *configList,
			       PDB *pdb, NamdState *state)
{
  this->molecule = molecule;
  this->parameters = parameters;
  this->simParameters = simParameters;
  this->configList = configList;
  this->pdb = pdb;
  this->state = state;
}

//======================================================================
// Private functions


#include "Node.bot.h"

/***************************************************************************
 * RCS INFORMATION:
 *
 *	$RCSfile: Node.C,v $
 *	$Author: jim $	$Locker:  $		$State: Exp $
 *	$Revision: 1.1002 $	$Date: 1997/02/13 04:43:10 $
 *
 ***************************************************************************
 * REVISION HISTORY:
 *
 * $Log: Node.C,v $
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
