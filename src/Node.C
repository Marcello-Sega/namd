/***************************************************************************/
/*                                                                         */
/*              (C) Copyright 1996 The Board of Trustees of the            */
/*                          University of Illinois                         */
/*                           All Rights Reserved                           */
/*                                                                         */
/***************************************************************************/
static char ident[] = "@(#)$Header: /home/cvs/namd/cvsroot/namd2/src/Node.C,v 1.4 1996/08/21 23:58:25 brunner Exp $";

#include <stdio.h>

#include "ckdefs.h"
#include "chare.h"
#include "c++interface.h"

#include "main.h"
#include "main.top.h"
#include "Node.top.h"
#include "Node.h"
#include "WorkDistrib.h"
#include "PatchMgr.h"
//#include "ProxyMgr.h"
//#include "MessageComm.h"
//#include "PatchMap.h"

//======================================================================
// Public Functions

//----------------------------------------------------------------------

//----------------------------------------------------------------------
// Node(void) can allocate some static data, or set all pointers to NULL

Node::Node(NodeInitMsg *msg)
{
  molecule = NULL;
  parameters = NULL;
  simParameters = NULL;
  configList = NULL;
  pdb = NULL;

  workDistrib = CLocalBranch(WorkDistrib,msg->workDistribGroup);
  patchMgr = CLocalBranch(PatchMgr,msg->patchMgrGroup);
}

//----------------------------------------------------------------------
// ~Node(void) needs to clean up everything.

Node::~Node(void)
{
}

//----------------------------------------------------------------------
// startup(void) receives all the startup data from the Namd object,
// and does everything necessary to run the simulation
//    Ask work_distrib for initial distribution
//    This must be split-phase, so results can be received from node 0.
void Node::startup(InitMsg *msg)
{

  CPrintf("Node %d startup\n",myid());

  workDistrib->parentNode(this);
  workDistrib->buildMapsFromScratch();
}


// run(void) runs the specified simulation for the specified number of
// steps, overriding the contents of the configuration file
void Node::run(void)
{
  CPrintf("Node::run() - invoked\n");
}

int Node::numNodes(void)
{
  return CNumPes();
}

int Node::myid(void)
{
  return CMyPe();
}

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
 *	$Author: brunner $	$Locker:  $		$State: Exp $
 *	$Revision: 1.4 $	$Date: 1996/08/21 23:58:25 $
 *
 ***************************************************************************
 * DESCRIPTION:
 *
 ***************************************************************************
 * REVISION HISTORY:
 *
 * $Log: Node.C,v $
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
 * Revision 1.7  1996/07/01 16:25:30  brunner
 * *** empty log message ***
 *
 * Revision 1.6  1996/06/14 18:39:00  brunner
 * *** empty log message ***
 *
 * Revision 1.5  1996/06/11 22:36:23  brunner
 * *** empty log message ***
 *
 * Revision 1.4  1996/06/10 22:04:14  brunner
 * *** empty log message ***
 *
 * Revision 1.3  1996/06/06 13:55:05  brunner
 * *** empty log message ***
 *
 * Revision 1.2  1996/06/04 16:12:18  brunner
 * Create dummy patches
 *
 * Revision 1.1  1996/05/30 20:16:09  brunner
 * Initial revision
 ***************************************************************************/
