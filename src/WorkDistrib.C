/***************************************************************************/
/*                                                                         */
/*              (C) Copyright 1996 The Board of Trustees of the            */
/*                          University of Illinois                         */
/*                           All Rights Reserved                           */
/*                                                                         */
/***************************************************************************/

/***************************************************************************/
/* DESCRIPTION:								   */
/*                                                                         */
/***************************************************************************/

static char ident[] = "@(#)$Header: /home/cvs/namd/cvsroot/namd2/src/WorkDistrib.C,v 1.7 1996/08/19 21:41:31 brunner Exp $";

#include <stdio.h>

#include "ckdefs.h"
#include "chare.h"
#include "c++interface.h"

#include "WorkDistrib.top.h"
#include "WorkDistrib.h"

#include "main.top.h"
#include "main.h"
#include "Node.h"
#include "PatchMap.h"
#include "ComputeMap.h"
#include "Vector.h"
#include "NamdTypes.h"

//======================================================================
// Public functions
//----------------------------------------------------------------------
WorkDistrib::WorkDistrib(InitMsg *msg)
{
}

//----------------------------------------------------------------------
WorkDistrib::~WorkDistrib(void)
{
  node = NULL;
}

//----------------------------------------------------------------------
void WorkDistrib::parentNode(Node *inode)
{
  node = inode;
}

//----------------------------------------------------------------------
void WorkDistrib::buildMapsFromScratch()
{
  int i;
  PatchMap *patchMap = &(node->patchMap);
  
  CPrintf("Building maps\n");
  mapPatches();
  mapComputes();

  node->patchMap.printPatchMap();
  //  node->computeMap.printComputeMap();

  // Send maps to other nodes.

  MapDistribMsg *mapMsg = new (MsgIndex(MapDistribMsg)) MapDistribMsg ;

  CBroadcastMsgBranch(WorkDistrib, saveMaps, mapMsg, thisgroup);

  // Create patches on this node

  for(i=0; i < patchMap->numPatches(); i++)
  {
    CPrintf("patchMgr->createPatch(%d,atoms,positions)\n",i);
  }

  // Move patches to the proper node
  for(i=0;i < patchMap->numPatches(); i++)
  {
    if (patchMap->node(i) != node->myid())
      CPrintf("patchMgr->movePatch(%d,%d)\n",i,patchMap->node(i));
  }
}

//----------------------------------------------------------------------
// saveMaps() is called when the map message is received
void WorkDistrib::saveMaps(MapDistribMsg *msg)
{
  if (node->myid() != 0)
  {
    CPrintf("Saving patch map, compute map\n");
  }
  else
  {
    CPrintf("Node 0 patch map already built\n");
  }
  CharmExit();
}


//======================================================================
// Private functions
//----------------------------------------------------------------------
void WorkDistrib::mapPatches(void)
{
  PatchMap *patchMap = &(node->patchMap);
  int xdim, ydim, zdim;
  int xi, yi, zi, pid;
  int i;
  Vector xmin, xmax;
  Position sysDim;

  CPrintf("Mapping patches\n");
  /*
   *  pdb->find_extremes(&xmin,&xmax);
   */

  xmax.x = xmax.y = xmax.z = 20.;
  xmin.x = xmin.y = xmin.z = 0.;

  sysDim.x = xmax.x - xmin.x;
  sysDim.y = xmax.y - xmin.y;
  sysDim.z = xmax.z - xmin.z;
  xdim = (int)((float)sysDim.x / patchSize);
  if ((xdim * patchSize) < sysDim.x)
    xdim++;

  ydim = (int)((float)sysDim.y / patchSize);
  if ((ydim * patchSize) < sysDim.y)
    ydim++;

  zdim = (int)((float)sysDim.z / patchSize);
  if ((zdim * patchSize) < sysDim.z)
    zdim++;

  patchMap->allocatePids(xdim, ydim, zdim);

  int assignedNode=0;
  for(i=0; i < patchMap->numPatches(); i++)
  {
    pid=patchMap->requestPid(&xi,&yi,&zi);
    patchMap->storePatch(pid,assignedNode,250, 
			 (xi*patchSize)/xdim+xmin.x,
			 (yi*patchSize)/ydim+xmin.y,
			 (zi*patchSize)/zdim+xmin.z,
			 (xi*patchSize)/xdim+xmin.x+patchSize,
			 (yi*patchSize)/ydim+xmin.y+patchSize,
			 (zi*patchSize)/zdim+xmin.z+patchSize);
    assignedNode++;
    if (node->numNodes()==assignedNode)
      assignedNode=0;
  }
}

//----------------------------------------------------------------------
void WorkDistrib::mapComputes(void)
{
  PatchMap *patchMap = &(node->patchMap);
  ComputeMap *computeMap = &(node->computeMap);

  CPrintf("Mapping computes\n");

  // We need to allocate computes for self, 1 and 2 away pairs for
  // electrostatics, and 1 angleForce for each node.  Then I might
  // throw in a few extras, in case I forget some.

  int numPotentialCids = 
    (patchMap->numPatches() * 125 + 1) / 2  + node->numNodes();

  computeMap->allocateCids(numPotentialCids);

  mapAngleComputes();

  mapElectComputes();
}

//----------------------------------------------------------------------
void WorkDistrib::mapAngleComputes()
{
  int i;
  PatchMap *patchMap = &(node->patchMap);
  ComputeMap *computeMap = &(node->computeMap);
  int numNodes = node->numNodes();

  // Figure out where to put the angleForce objects, and hook up the
  // dependencies.

  ComputeID *cid = new ComputeID[numNodes];

  for(i=0; i<node->numNodes(); i++)
  {
    cid[i]=computeMap->storeCompute(i,patchMap->numPatches(),
				    angleForceType);
  }

  PatchID j;

  for(j=0;j<patchMap->numPatches();j++)
  {
    patchMap->newCid(j,cid[patchMap->node(j)]);
    computeMap->newPid(cid[patchMap->node(j)],j);
  }

  delete cid;
}

//----------------------------------------------------------------------
void WorkDistrib::mapElectComputes(void)
{
  // For each patch, create 1 electrostatic object for self-interaction.
  // Then create 1 for each 1-away and 2-away neighbor which has a larger
  // pid.

  PatchMap *patchMap = &(node->patchMap);
  ComputeMap *computeMap = &(node->computeMap);

  PatchID oneAway[PatchMap::MaxOneAway];
  PatchID twoAway[PatchMap::MaxTwoAway];

  PatchID i;
  ComputeID cid;
  int numNeighbors;
  int j;

  for(i=0; i<patchMap->numPatches(); i++)
  {
    // self-interaction
    cid=computeMap->storeCompute(patchMap->node(i),1,electForceType);
    computeMap->newPid(cid,i);
    patchMap->newCid(i,cid);

    // one-away neighbors
    numNeighbors=patchMap->oneAwayNeighbors(i,oneAway);
    for(j=0;j<numNeighbors;j++)
    {
      if (i < oneAway[j])
      {
	cid=computeMap->storeCompute(patchMap->node(i),2,
				     electForceType);
	computeMap->newPid(cid,i);
	computeMap->newPid(cid,oneAway[j]);
	patchMap->newCid(i,cid);
	patchMap->newCid(oneAway[j],cid);
      }
    }

    // two-away neighbors
    numNeighbors=patchMap->twoAwayNeighbors(i,oneAway);
    for(j=0;j<numNeighbors;j++)
    {
      if (i < oneAway[j])
      {
	cid=computeMap->storeCompute(patchMap->node(i),2,
				     electForceType);
	computeMap->newPid(cid,i);
	computeMap->newPid(cid,oneAway[j]);
	patchMap->newCid(i,cid);
	patchMap->newCid(oneAway[j],cid);
      }
    }
  }
}

#include "WorkDistrib.bot.h"

/***************************************************************************
 * RCS INFORMATION:
 *
 *	$RCSfile: WorkDistrib.C,v $
 *	$Author: brunner $	$Locker:  $		$State: Exp $
 *	$Revision: 1.7 $	$Date: 1996/08/19 21:41:31 $
 *
 ***************************************************************************
 * REVISION HISTORY:
 *
 * $Log: WorkDistrib.C,v $
 * Revision 1.7  1996/08/19 21:41:31  brunner
 * *** empty log message ***
 *
 * Revision 1.6  1996/08/19 21:37:02  brunner
 * Create Patches from PDB data
 *
 * Revision 1.5  1996/08/19 20:39:11  brunner
 * *** empty log message ***
 *
 * Revision 1.4  1996/08/16 21:55:56  brunner
 * *** empty log message ***
 *
 * Revision 1.3  1996/08/16 21:41:11  brunner
 * *** empty log message ***
 *
 * Revision 1.2  1996/08/16 21:16:04  ari
 * *** empty log message ***
 *
 * Revision 1.1  1996/08/16 20:34:23  brunner
 * Initial revision
 *
 * Revision 1.6  1996/08/03 20:08:09  brunner
 * *** empty log message ***
 *
 * Revision 1.5  1996/07/01 16:25:30  brunner
 * *** empty log message ***
 *
 * Revision 1.4  1996/06/14 18:39:00  brunner
 * *** empty log message ***
 *
 * Revision 1.3  1996/06/11 22:36:23  brunner
 * *** empty log message ***
 *
 * Revision 1.2  1996/06/04 16:12:18  brunner
 * Create dummy patches
 *
 * Revision 1.1  1996/05/30 20:16:09  brunner
 * Initial revision
 *
 * Revision 1.1  1996/05/30 20:11:09  brunner
 * Initial revision
 *
 ***************************************************************************/

