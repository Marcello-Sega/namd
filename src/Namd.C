/***************************************************************************/
/*                                                                         */
/*              (C) Copyright 1996 The Board of Trustees of the            */
/*                          University of Illinois                         */
/*                           All Rights Reserved                           */
/*                                                                         */
/***************************************************************************/

static char ident[] = "@(#)$Header: /home/cvs/namd/cvsroot/namd2/src/Attic/Namd.C,v 1.7 1996/08/29 00:50:42 ari Exp $";

#include "ckdefs.h"
#include "chare.h"
#include "c++interface.h"

#include "main.top.h"
#include "main.h"
#include "Namd.h"
#include "Molecule.h"
#include "Parameters.h"
#include "SimParameters.h"
#include "ConfigList.h"
#include "Node.top.h"
#include "Node.h"
#include "WorkDistrib.top.h"
#include "WorkDistrib.h"
#include "PatchMgr.top.h"
#include "PatchMgr.h"

// Namd(void ) is the constructor for the startup node.  It needs to
// read in file data,
Namd::Namd(void)
{
  InitMsg *initmsg;
  PatchMgrInitMsg *pminitmsg;

  NodeInitMsg *node_msg = new (MsgIndex(NodeInitMsg)) NodeInitMsg;
  
  // Create WorkDistrib and send it an empty message
  initmsg = new (MsgIndex(InitMsg)) InitMsg;
  node_msg->workDistribGroup
    = new_group(WorkDistrib, initmsg);

  // Create PatchMgr and send it an empty message
  pminitmsg = new (MsgIndex(PatchMgrInitMsg)) PatchMgrInitMsg;
  pminitmsg->workDistribGroup = node_msg->workDistribGroup
    = new_group(PatchMgr, pminitmsg);

  // Create the Node object and send it the IDs of all the other
  // parallel objects.
  nodeGroup = new_group(Node, node_msg);
  node = CLocalBranch(Node,nodeGroup);
}


// ~Namd(void) just needs to tell all the slave nodes to die.
Namd::~Namd(void)
{
}


// startup(char *) 
void Namd::startup(char *confFile)
{
  InitMsg *initmsg;

  namdState.configFileInit(confFile);
  if (namdState.status()) {
    CPrintf("Namd::startup() - could not initialize namdState from %s\n", 
      confFile);
    CharmExit();
  }

  // Give node[0] pointers to the data objects, so it can use them,
  // or send them on as messages elsewhere.
  node->saveMolDataPointers(namdState.molecule,namdState.parameters,
			    namdState.simParameters,namdState.configList,
			    namdState.pdb);

  // Tell Node to do any startup work
  initmsg = new (MsgIndex(InitMsg)) InitMsg;
  CSendMsgBranch(Node, startup, initmsg, nodeGroup, CMyPe());

}


// run(void) runs the specified simulation to completion
void Namd::run(void)
{
  CPrintf("Namd::run() invoked\n");
//  node->run();
  CharmExit();
}

/***************************************************************************
 * RCS INFORMATION:
 *
 *	$RCSfile: Namd.C,v $
 *	$Author: ari $	$Locker:  $		$State: Exp $
 *	$Revision: 1.7 $	$Date: 1996/08/29 00:50:42 $
 *
 ***************************************************************************
 * DESCRIPTION:
 *
 ***************************************************************************
 * REVISION HISTORY:
 *
 * $Log: Namd.C,v $
 * Revision 1.7  1996/08/29 00:50:42  ari
 * *** empty log message ***
 *
 * Revision 1.6  1996/08/21 23:58:25  brunner
 * *** empty log message ***
 *
 * Revision 1.5  1996/08/19 22:05:31  ari
 * *** empty log message ***
 *
 * Revision 1.4  1996/08/16 21:19:34  ari
 * *** empty log message ***
 *
 * Revision 1.3  1996/08/16 20:33:09  brunner
 * Not Working
 *
 * Revision 1.2  1996/08/16 04:39:46  ari
 * *** empty log message ***
 *
 * Revision 1.1  1996/08/16 01:54:59  ari
 * Initial revision
 *
 ***************************************************************************/
