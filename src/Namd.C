/***************************************************************************/
/*                                                                         */
/*              (C) Copyright 1996 The Board of Trustees of the            */
/*                          University of Illinois                         */
/*                           All Rights Reserved                           */
/*                                                                         */
/***************************************************************************/

static char ident[] = "@(#)$Header: /home/cvs/namd/cvsroot/namd2/src/Attic/Namd.C,v 1.12 1996/12/05 21:37:43 ari Exp $";

#include "unistd.h"

#include "ckdefs.h"
#include "chare.h"
#include "c++interface.h"

#include "main.top.h"
#include "main.h"
#include "BOCgroup.h"
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
#include "ComputeMgr.top.h"
#include "ComputeMgr.h"
#include "ProxyMgr.top.h"
#include "ProxyMgr.h"


// Namd(void ) is the constructor for the startup node.  It needs to
// read in file data,
Namd::Namd(void)
{
  InitMsg *initmsg;
  BOCgroup group;

  // Create WorkDistrib and send it an empty message
  initmsg = new (MsgIndex(InitMsg)) InitMsg;
  group.workDistrib = new_group(WorkDistrib, initmsg);

  // Create ProxyMgr
  initmsg = new (MsgIndex(InitMsg)) InitMsg;
  group.proxyMgr = new_group(ProxyMgr, initmsg);

  // Create PatchMgr
  initmsg = new (MsgIndex(InitMsg)) InitMsg;
  group.patchMgr = new_group(PatchMgr, initmsg);

  // Create ComputeMgr
  initmsg = new (MsgIndex(InitMsg)) InitMsg;
  group.computeMgr = new_group(ComputeMgr, initmsg);

  // Create the Node object and send it the IDs of all the other
  // parallel objects.
  GroupInitMsg *msg = new (MsgIndex(GroupInitMsg)) GroupInitMsg;
  msg->group = group;

  nodeGroup = new_group(Node, msg);
}


// ~Namd(void) just needs to tell all the slave nodes to die.
Namd::~Namd(void)
{
}


// startup(char *) 
void Namd::startup(char *confFile)
{
  namdState.configFileInit(confFile);
  if (namdState.status()) {
    CPrintf("Namd::startup() - could not initialize namdState from %s\n", 
      confFile);
    CharmExit();
  }

  // Give node[0] pointers to the data objects, so it can use them,
  // or send them on as messages elsewhere.
  Node::Object()->saveMolDataPointers(namdState.molecule,namdState.parameters,
			    namdState.simParameters,namdState.configList,
			    namdState.pdb);

  Node::messageStartup();
}



/***************************************************************************
 * RCS INFORMATION:
 *
 *	$RCSfile: Namd.C,v $
 *	$Author: ari $	$Locker:  $		$State: Exp $
 *	$Revision: 1.12 $	$Date: 1996/12/05 21:37:43 $
 *
 ***************************************************************************
 * DESCRIPTION:
 *
 ***************************************************************************
 * REVISION HISTORY:
 *
 * $Log: Namd.C,v $
 * Revision 1.12  1996/12/05 21:37:43  ari
 * *** empty log message ***
 *
 * Revision 1.11  1996/11/30 00:38:14  jim
 * added ComputeMgr creation
 *
 * Revision 1.10  1996/11/22 00:18:51  ari
 * *** empty log message ***
 *
 * Revision 1.9  1996/11/05 16:59:58  ari
 * *** empty log message ***
 *
 * Revision 1.8  1996/09/03 22:54:09  ari
 * *** empty log message ***
 *
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
