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

static char ident[] = "@(#)$Header: /home/cvs/namd/cvsroot/namd2/src/Node.C,v 1.16 1996/11/30 00:44:24 jim Exp $";


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
#include "AtomMap.h"
#include "Sequencer.h"
//#include "ProxyMgr.h"
//#include "MessageComm.h"
//#include "PatchMap.h"

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

  molecule = NULL;
  parameters = NULL;
  simParameters = NULL;
  configList = NULL;
  pdb = NULL;
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
}

//----------------------------------------------------------------------
// ~Node(void) needs to clean up everything.

Node::~Node(void)
{
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
  CSendMsgBranch(Node, startup, msg, group.node, CMyPe() );
}

void Node::startup(InitMsg *msg)
{
  AtomMap::Object()->allocateMap(molecule->numAtoms);

  delete msg;

  DebugM(1, "Node::startup() Pe=" << CMyPe() << "\n");

  // Thread initialization
  if (CthImplemented()) {
    CthSetStrategyDefault(CthSelf());
  } else {
    DebugM(1, "Node::startup() - Oh no, tiny elvis, threads not implemented\n");
    CharmExit();
  }

  workDistrib = CLocalBranch(WorkDistrib,group.workDistrib);

  workDistrib->buildMaps();
  workDistrib->sendMaps();
  workDistrib->awaitMaps();

  ComputeMap::Object()->printComputeMap();

  workDistrib->createPatches();
  // workDistrib->createComputes();

  computeMgr = CLocalBranch(ComputeMgr,group.computeMgr);
  DebugM(3, "Trying to create computes.\n");
  computeMgr->createComputes(ComputeMap::Object());

  patchMgr = CLocalBranch(PatchMgr,group.patchMgr);

  messageStartupDone();   // collect on master node
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

  HomePatchList *hpl = PatchMap::Object()->homePatchList();
  ResizeArrayIter<HomePatchElem> ai(*hpl);

  DebugM(4, "Node::run() - iterating over home patches!\n");
  for (ai=ai.begin(); ai != ai.end(); ai++) {
    HomePatch *p = (*ai).p;
    DebugM(1, "Node::run() - signaling patch "<< p->getPatchID() << endl);
    Sequencer *sequencer = new Sequencer(p);
    p->useSequencer(sequencer);
    p->runSequencer(1);
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
			       PDB *pdb)
{
  this->molecule = molecule;
  this->parameters = parameters;
  this->simParameters = simParameters;
  this->configList = configList;
  this->pdb = pdb;
}

//======================================================================
// Private functions


#include "Node.bot.h"

/***************************************************************************
 * RCS INFORMATION:
 *
 *	$RCSfile: Node.C,v $
 *	$Author: jim $	$Locker:  $		$State: Exp $
 *	$Revision: 1.16 $	$Date: 1996/11/30 00:44:24 $
 *
 ***************************************************************************
 * REVISION HISTORY:
 *
 * $Log: Node.C,v $
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
