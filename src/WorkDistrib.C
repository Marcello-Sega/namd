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

static char ident[] = "@(#)$Header: /home/cvs/namd/cvsroot/namd2/src/WorkDistrib.C,v 1.779 1997/02/06 15:53:32 ari Exp $";

#include <stdio.h>

#include "ckdefs.h"
#include "chare.h"
#include "c++interface.h"

#include "BOCgroup.h"
#include "WorkDistrib.top.h"
#include "WorkDistrib.h"

#include "main.top.h"
#include "main.h"
#include "Node.h"
#include "PatchMgr.h"
#include "PatchMap.h"
#include "ComputeMap.h"
#include "Compute.h"
#include "Vector.h"
#include "NamdTypes.h"
#include "PDB.h"
#include "SimParameters.h"

#define MIN_DEBUG_LEVEL 3
#define DEBUGM
#include "Debug.h"

//======================================================================
// Public functions
//----------------------------------------------------------------------
WorkDistrib::WorkDistrib(InitMsg *msg)
{
  DebugM(4,"WorkDistrib::WorkDistrib() - constructing\n");
  mapsArrived = false;
  awaitingMaps = false;
  DebugM(4,"WorkDistrib::WorkDistrib() - done constructing\n");
}

//----------------------------------------------------------------------
WorkDistrib::~WorkDistrib(void)
{
}

//----------------------------------------------------------------------
void WorkDistrib::buildMaps(void)
{
  int i;

  DebugM(4,"Building maps\n");
  mapPatches();
  mapComputes();
}

//----------------------------------------------------------------------
void WorkDistrib::sendMaps(void)
{
  MapDistribMsg *mapMsg = new (MsgIndex(MapDistribMsg)) MapDistribMsg ;

  mapMsg->patchMap = PatchMap::Object();
  mapMsg->computeMap = ComputeMap::Object();

  DebugM(3,"calling CBroadcastMsgBranch\n");
  CBroadcastMsgBranch(WorkDistrib, saveMaps, mapMsg, thisgroup);
  DebugM(3,"called CBroadcastMsgBranch\n");
  mapsArrived = true;
}

//----------------------------------------------------------------------
void WorkDistrib::createComputes(void)
{
  DebugM(7,"I don't know how to create computes yet\n");
  CharmExit();
}

//----------------------------------------------------------------------
void WorkDistrib::createPatches(void)
{
  int i,j;

  PatchMap *patchMap = PatchMap::Object();
  Node *node = CLocalBranch(Node,group.node);
  PatchMgr *patchMgr = CLocalBranch(PatchMgr,group.patchMgr);
  SimParameters *params = node->simParameters;
  PDB *pdb = node->pdb;

  int numPatches = patchMap->numPatches();
  int numAtoms = pdb->num_atoms();

  Vector *positions = new Position[numAtoms];
  pdb->get_all_positions(positions);

  AtomIDList *atomIDs = new AtomIDList[numPatches];
  PositionList *atomPositions = new PositionList[numPatches];
  VelocityList *atomVelocities = new VelocityList[numPatches];

  Lattice lattice = params->lattice;
  Velocity vel(0.,0.,0.);

  for(i=0; i < numAtoms; i++)
  {
    if ( ! ( i % 1000 ) )
    {
      DebugM(4,"Assigned " << i << " atoms to patches so far.\n");
    }
    int pid = patchMap->assignToPatch(positions[i]);
    atomIDs[pid].add(i);
    atomPositions[pid].add(positions[i]);
    atomVelocities[pid].add(vel);
  }

  delete [] positions;

  for(i=0; i < numPatches; i++)
  {
    if ( ! ( i % 100 ) )
    {
      DebugM(4,"Created " << i << " patches so far.\n");
    }

    Position center(0.5*(patchMap->minX(i)+patchMap->maxX(i)),
		    0.5*(patchMap->minY(i)+patchMap->maxY(i)),
		    0.5*(patchMap->minZ(i)+patchMap->maxZ(i)));

    for(int j=0; j < atomIDs[i].size(); j++)
    {
      atomPositions[i][j] = (lattice.nearest(atomPositions[i][j],center));
    }

    patchMgr->createHomePatch(i,atomIDs[i],atomPositions[i],atomVelocities[i]);
  }

  delete [] atomIDs;
  delete [] atomPositions;
  delete [] atomVelocities;

  // Move patches to the proper node
  for(i=0;i < patchMap->numPatches(); i++)
  {
    if (patchMap->node(i) != node->myid() )
    {
      DebugM(3,"patchMgr->movePatch("
	<< i << "," << patchMap->node(i) << ")\n");
      patchMgr->movePatch(i,patchMap->node(i));
    }
  }
  patchMgr->sendMovePatches();
}

//----------------------------------------------------------------------
// saveMaps() is called when the map message is received
void WorkDistrib::saveMaps(MapDistribMsg *msg)
{
  Node *node = CLocalBranch(Node,group.node);

  if (node->myid() != 0)
  {
    DebugM(4,"Saving patch map, compute map\n");
  }
  else
  {
    DebugM(4,"Node 0 patch map built\n");
  }

  mapsArrived = true;
}

//----------------------------------------------------------------------
// awaitMaps() is called when node needs to wait for the map message
void WorkDistrib::awaitMaps()
{
  while (!mapsArrived)
  {
    awaitingMapsTh = CthSelf();
    awaitingMaps = true;
    DebugM(4,"suspending in awaitMaps(), thread = " << awaitingMapsTh << "\n");
    CthYield();
  }
}


//======================================================================
// Private functions
//----------------------------------------------------------------------
void WorkDistrib::mapPatches(void)
{
  PatchMap *patchMap = PatchMap::Object();
  Node *node = CLocalBranch(Node, group.node);
  SimParameters *params = node->simParameters;

  int xdim, ydim, zdim;
  int xper, yper, zper;
  int xi, yi, zi, pid;
  int i;
  Vector xmin, xmax;
  Vector sysDim;

  DebugM(4,"Mapping patches\n");
  node->pdb->find_extremes(&xmin,&xmax);

  DebugM(4,"xmin.x = " << xmin.x << endl);
  DebugM(4,"xmin.y = " << xmin.y << endl);
  DebugM(4,"xmin.z = " << xmin.z << endl);

  DebugM(4,"xmax.x = " << xmax.x << endl);
  DebugM(4,"xmax.y = " << xmax.y << endl);
  DebugM(4,"xmax.z = " << xmax.z << endl);

  if ( params->cellBasisVector1.x )
  {
    xper = 1;
    sysDim.x = params->cellBasisVector1.x;
    xdim = (int)(sysDim.x / patchSize);
    if ( xdim < 2 ) xdim = 2;
  }
  else
  {
    xper = 0;
    sysDim.x = xmax.x - xmin.x;
    xdim = (int)((float)sysDim.x / patchSize);
    if ((xdim * patchSize) < sysDim.x)
      xdim++;
    sysDim.x = xdim * patchSize;
  }

  if ( params->cellBasisVector2.y )
  {
    yper = 1;
    sysDim.y = params->cellBasisVector2.y;
    ydim = (int)(sysDim.y / patchSize);
    if ( ydim < 2 ) ydim = 2;
  }
  else
  {
    yper = 0;
    sysDim.y = xmax.y - xmin.y;
    ydim = (int)((float)sysDim.y / patchSize);
    if ((ydim * patchSize) < sysDim.y)
      ydim++;
    sysDim.y = ydim * patchSize;
  }

  if ( params->cellBasisVector3.z )
  {
    zper = 1;
    sysDim.z = params->cellBasisVector3.z;
    zdim = (sysDim.z / patchSize);
    if ( zdim < 2 ) zdim = 2;
  }
  else
  {
    zper = 0;
    sysDim.z = xmax.z - xmin.z;
    zdim = (int)((float)sysDim.z / patchSize);
    if ((zdim * patchSize) < sysDim.z)
      zdim++;
    sysDim.z = zdim * patchSize;
  }

  patchMap->setPeriodicity(xper,yper,zper);
  patchMap->allocatePids(xdim, ydim, zdim);
  patchMap->setGridOriginAndLength(xmin,sysDim);

  int assignedNode=0;
  for(i=0; i < patchMap->numPatches(); i++)
  {
    pid=patchMap->requestPid(&xi,&yi,&zi);
    patchMap->storePatch(pid,assignedNode,250, 
			 ((float)xi/(float)xdim)*sysDim.x+xmin.x,
			 ((float)yi/(float)ydim)*sysDim.y+xmin.y,
			 ((float)zi/(float)zdim)*sysDim.z+xmin.z,
			 ((float)(xi+1)/(float)xdim)*sysDim.x+xmin.x,
			 ((float)(yi+1)/(float)ydim)*sysDim.y+xmin.y,
			 ((float)(zi+1)/(float)zdim)*sysDim.z+xmin.z);
    assignedNode++;
    if (node->numNodes()==assignedNode)
      assignedNode=0;
  }
}

//----------------------------------------------------------------------
void WorkDistrib::mapComputes(void)
{
  PatchMap *patchMap = PatchMap::Object();
  ComputeMap *computeMap = ComputeMap::Object();
  Node *node = CLocalBranch(Node, group.node);

  DebugM(4,"Mapping computes\n");

  // We need to allocate computes for self, 1 and 2 away pairs for
  // electrostatics, and 1 angleForce for each node.  Then I might
  // throw in a few extras, in case I forget some.

  int numPotentialCids = 
    patchMap->numPatches() * (124/2+1) + node->numNodes() * 10;

  computeMap->allocateCids(numPotentialCids);

  // Handle full electrostatics
  if ( node->simParameters->fullDirectOn )
    mapComputeHomePatches(computeFullDirectType);
  if ( node->simParameters->FMAOn )
    mapComputeHomePatches(computeDPMTAType);

  mapComputeNonbonded();
  mapComputeHomePatches(computeNonbondedExclType);
  mapComputeHomePatches(computeBondsType);
  mapComputeHomePatches(computeAnglesType);
  mapComputeHomePatches(computeDihedralsType);
  mapComputeHomePatches(computeImpropersType);
}

//----------------------------------------------------------------------
void WorkDistrib::mapComputeHomePatches(ComputeType type)
{
  PatchMap *patchMap = PatchMap::Object();
  ComputeMap *computeMap = ComputeMap::Object();
  Node *node = CLocalBranch(Node, group.node);

  int i;

  int numNodes = node->numNodes();

  ComputeID *cid = new ComputeID[numNodes];

  for(i=0; i<node->numNodes(); i++)
  {
    cid[i]=computeMap->storeCompute(i,patchMap->numPatches(),type);
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
void WorkDistrib::mapComputeNonbonded(void)
{
  // For each patch, create 1 electrostatic object for self-interaction.
  // Then create 1 for each 1-away and 2-away neighbor which has a larger
  // pid.

  PatchMap *patchMap = PatchMap::Object();
  ComputeMap *computeMap = ComputeMap::Object();

  PatchID oneAway[PatchMap::MaxOneAway];
  PatchID oneAwayTrans[PatchMap::MaxOneAway];
  PatchID twoAway[PatchMap::MaxTwoAway];

  PatchID i;
  ComputeID cid;
  int numNeighbors;
  int j;

  for(i=0; i<patchMap->numPatches(); i++)
  {
    // self-interaction
    cid=computeMap->storeCompute(patchMap->node(i),1,computeNonbondedSelfType);
    computeMap->newPid(cid,i);
    patchMap->newCid(i,cid);

    // one-away neighbors
    numNeighbors=patchMap->oneAwayNeighbors(i,oneAway,oneAwayTrans);
    for(j=0;j<numNeighbors;j++)
    {
      if (i < oneAway[j])
      {
	cid=computeMap->storeCompute(patchMap->node(i),2,
				     computeNonbondedPairType);
	computeMap->newPid(cid,i);
	computeMap->newPid(cid,oneAway[j],oneAwayTrans[j]);
	patchMap->newCid(i,cid);
	patchMap->newCid(oneAway[j],cid);
      }
    }

/*
    // two-away neighbors
    numNeighbors=patchMap->twoAwayNeighbors(i,twoAway);
    for(j=0;j<numNeighbors;j++)
    {
      if (i < twoAway[j])
      {
	cid=computeMap->storeCompute(patchMap->node(i),2,
				     computeNonbondedPairType);
	computeMap->newPid(cid,i);
	computeMap->newPid(cid,twoAway[j]);
	patchMap->newCid(i,cid);
	patchMap->newCid(twoAway[j],cid);
      }
    }
*/
  }
}

//----------------------------------------------------------------------
void WorkDistrib::messageEnqueueWork(Compute *compute) {
  DebugM(2,"WorkDistrib::messageEnqueueWork() triggered\n");
  LocalWorkMsg *msg = new (MsgIndex(LocalWorkMsg)) LocalWorkMsg;
  msg->compute = compute; // pointer is valid since send is to local Pe
  CSendMsgBranch(WorkDistrib, enqueueWork, msg, group.workDistrib, CMyPe() );
}

void WorkDistrib::enqueueWork(LocalWorkMsg *msg) {
  msg->compute->doWork();
  delete msg;
}

void WorkDistrib::messageMovePatchDone() {
  // Send msg to WorkDistrib that all patchMoves are completed
  DoneMsg *msg = new (MsgIndex(DoneMsg)) DoneMsg;
  CSendMsgBranch(WorkDistrib, movePatchDone, msg, group.workDistrib, CMyPe());
}

void WorkDistrib::movePatchDone(DoneMsg *msg) {
  delete msg;
}

#include "WorkDistrib.bot.h"

/***************************************************************************
 * RCS INFORMATION:
 *
 *	$RCSfile: WorkDistrib.C,v $
 *	$Author: ari $	$Locker:  $		$State: Exp $
 *	$Revision: 1.779 $	$Date: 1997/02/06 15:53:32 $
 *
 ***************************************************************************
 * REVISION HISTORY:
 *
 * $Log: WorkDistrib.C,v $
 * Revision 1.779  1997/02/06 15:53:32  ari
 * Updating Revision Line, getting rid of branches
 *
 * Revision 1.778.2.2  1997/02/06 05:00:45  jim
 * Added creation of full electrostatics objects.
 *
 * Revision 1.778.2.1  1997/02/06 02:35:35  jim
 * Implemented periodic boundary conditions - may not work with
 * atom migration yet, but doesn't seem to alter calculation,
 * appears to work correctly when turned on.
 * NamdState chdir's to same directory as config file in argument.
 *
 * Revision 1.778  1997/01/28 00:31:29  ari
 * internal release uplevel to 1.778
 *
 * Revision 1.777.2.1  1997/01/22 21:42:13  jim
 * Larger patches, no two-away computes, small tweak to inner loop.
 *
 * Revision 1.777  1997/01/17 19:37:04  ari
 * Internal CVS leveling release.  Start development code work
 * at 1.777.1.1.
 *
 * Revision 1.28  1996/12/27 22:22:33  nealk
 * Added debugging code.
 *
 * Revision 1.27  1996/12/19 02:26:30  jim
 * Node::startup2 is now triggered by quiescence
 *
 * Revision 1.26  1996/12/13 08:54:10  jim
 * now moves patches
 *
 * Revision 1.25  1996/12/12 20:14:50  milind
 * *** empty log message ***
 *
 * Revision 1.24  1996/12/12 17:24:04  jim
 * forgot to fill in message fields
 *
 * Revision 1.23  1996/12/12 08:58:58  jim
 * added movePatch call, changed CthResume to CthAwaken
 *
 * Revision 1.22  1996/12/01 21:02:37  jim
 * now adds all existing compute objects to map
 *
 * Revision 1.21  1996/12/01 02:29:59  jim
 * improved number of computes estimate
 *
 * Revision 1.20  1996/11/30 21:59:12  jim
 * fixed use of wrong array in nonbonded allocation
 *
 * Revision 1.19  1996/11/30 21:09:15  jim
 * cleaned up debug messages
 *
 * Revision 1.18  1996/11/30 01:27:34  jim
 * switched to realistic ComputeType definitions
 *
 * Revision 1.17  1996/11/22 01:02:18  ari
 * *** empty log message ***
 *
 * Revision 1.16  1996/11/22 00:18:51  ari
 * *** empty log message ***
 *
 * Revision 1.15  1996/11/05 16:59:58  ari
 * *** empty log message ***
 *
 * Revision 1.14  1996/10/29 23:35:27  ari
 * *** empty log message ***
 *
 * Revision 1.13  1996/10/29 17:18:20  brunner
 * Added initial thread stuff, but it only works on 1 processor
 *
 * Revision 1.12  1996/10/16 08:22:39  ari
 * *** empty log message ***
 *
 * Revision 1.11  1996/10/04 22:23:34  brunner
 * Empty method for createComputes
 *
 * Revision 1.10  1996/09/03 22:45:35  ari
 * *** empty log message ***
 *
 * Revision 1.9  1996/08/23 22:03:52  brunner
 * *** empty log message ***
 *
 * Revision 1.8  1996/08/21 23:58:25  brunner
 * *** empty log message ***
 *
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

