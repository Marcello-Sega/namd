/***************************************************************************/
/*                                                                         */
/*              (C) Copyright 1996 The Board of Trustees of the            */
/*                          University of Illinois                         */
/*                           All Rights Reserved                           */
/*                                                                         */
/***************************************************************************/

static char ident[] = "@(#)$Header: /home/cvs/namd/cvsroot/namd2/src/Attic/Namd.C,v 1.778 1997/01/28 00:30:55 ari Exp $";

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
#include "ReductionMgr.top.h"
#include "ReductionMgr.h"
#include "CollectionMgr.top.h"
#include "CollectionMgr.h"
#include "CollectionMaster.top.h"
#include "CollectionMaster.h"


// Namd(void ) is the constructor for the startup node.  It needs to
// read in file data,
Namd::Namd(void)
{
  BOCgroup group;

  // Create WorkDistrib and send it an empty message
  InitMsg *initmsg1 = new (MsgIndex(InitMsg)) InitMsg;
  group.workDistrib = new_group(WorkDistrib, initmsg1);

  // Create ProxyMgr
  InitMsg *initmsg2 = new (MsgIndex(InitMsg)) InitMsg;
  group.proxyMgr = new_group(ProxyMgr, initmsg2);

  // Create PatchMgr
  InitMsg *initmsg3 = new (MsgIndex(InitMsg)) InitMsg;
  group.patchMgr = new_group(PatchMgr, initmsg3);

  // Create ComputeMgr
  InitMsg *initmsg4 = new (MsgIndex(InitMsg)) InitMsg;
  group.computeMgr = new_group(ComputeMgr, initmsg4);

  // Create ReductionMgr
  InitMsg *initmsg5 = new (MsgIndex(InitMsg)) InitMsg;
  group.reductionMgr = new_group(ReductionMgr, initmsg5);

  // Create Collection system
  InitMsg *initmsg6 = new (MsgIndex(InitMsg)) InitMsg;
  ChareIDType collectionMaster;
  new_chare2(CollectionMaster,initmsg6,&collectionMaster,0);
  SlaveInitMsg *initmsg7 = new (MsgIndex(SlaveInitMsg)) SlaveInitMsg;
  initmsg7->master = collectionMaster;
  group.collectionMgr = new_group(CollectionMgr,initmsg7);

  // Create the Node object and send it the IDs of all the other
  // parallel objects.
  GroupInitMsg *msg = new (MsgIndex(GroupInitMsg)) GroupInitMsg;
  msg->group = group;

  nodeGroup = new_group(Node, msg);
}


// ~Namd(void) just needs to tell all the slave nodes to die.
Namd::~Namd(void)
{
  CPrintf("Namd::~Namd() called\n");
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
			    namdState.pdb,&namdState);

  Node::messageStartup();
}



/***************************************************************************
 * RCS INFORMATION:
 *
 *	$RCSfile: Namd.C,v $
 *	$Author: ari $	$Locker:  $		$State: Exp $
 *	$Revision: 1.778 $	$Date: 1997/01/28 00:30:55 $
 *
 ***************************************************************************
 * DESCRIPTION:
 *
 ***************************************************************************
 * REVISION HISTORY:
 *
 * $Log: Namd.C,v $
 * Revision 1.778  1997/01/28 00:30:55  ari
 * internal release uplevel to 1.778
 *
 * Revision 1.777  1997/01/17 19:36:30  ari
 * Internal CVS leveling release.  Start development code work
 * at 1.777.1.1.
 *
 * Revision 1.16  1997/01/13 23:26:19  jim
 * create collection system
 *
 * Revision 1.15  1997/01/09 20:48:10  jim
 * added Controller code
 *
 * Revision 1.14  1996/12/26 22:26:29  nealk
 * Added instantiation of reduction manager.
 *
 * Revision 1.13  1996/12/10 00:13:12  ari
 * *** empty log message ***
 *
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
