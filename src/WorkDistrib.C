/***************************************************************************/
/*         (C) Copyright 1996,1997 The Board of Trustees of the            */
/*                          University of Illinois                         */
/*                           All Rights Reserved                           */
/***************************************************************************/

/***************************************************************************
 * DESCRIPTION:	Currently, WorkDistrib generates the layout of the Patches,
 *              directs the construction and distribution of Computes and
 *              associates Computes with Patches.
 *                                                                         
 ***************************************************************************/

static char ident[] = "@(#)$Header: /home/cvs/namd/cvsroot/namd2/src/WorkDistrib.C,v 1.1035 1997/09/30 16:22:19 brunner Exp $";

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
#include "NamdTypes.h"
#include "PDB.h"
#include "SimParameters.h"
#include "Molecule.h"
#include "NamdOneTools.h"
#include "Compute.h"
#include "ComputeMap.h"
#include "RecBisection.h"

#define MIN_DEBUG_LEVEL 4
#define DEBUGM
#include "Debug.h"

extern "C" long int lrand48(void);

//======================================================================
// Public functions
//----------------------------------------------------------------------
WorkDistrib::WorkDistrib(InitMsg *msg)
{
  delete msg;

  group.workDistrib = thisgroup;
  mapsArrived = false;
  awaitingMaps = false;
}

//----------------------------------------------------------------------
WorkDistrib::~WorkDistrib(void)
{ }


//----------------------------------------------------------------------
void WorkDistrib::sendMaps(void)
{
  MapDistribMsg *mapMsg = new (MsgIndex(MapDistribMsg)) MapDistribMsg ;

  mapMsg->patchMap = PatchMap::Object();
  mapMsg->computeMap = ComputeMap::Object();

  CBroadcastMsgBranch(WorkDistrib, saveMaps, mapMsg, thisgroup);
  mapsArrived = true;
}

//----------------------------------------------------------------------
void WorkDistrib::saveComputeMap(int ep, int chareID)
{
  saveComputeMapReturnEP = ep;
  saveComputeMapReturnChareID = chareID;
  saveComputeMapCount = CNumPes();

  ComputeMap *computeMap = ComputeMap::Object();

  for (int i=0; i<computeMap->numComputes(); i++) {
    DebugM(3, "ComputeMap (" << i << ") node = " << computeMap->node(i) << " newNode = " << computeMap->newNode(i) << "\n");
  }
  
  ComputeMapDistribMsg *mapMsg 
    = new (MsgIndex(ComputeMapDistribMsg)) ComputeMapDistribMsg ;

  mapMsg->computeMap = ComputeMap::Object();

  CBroadcastMsgBranch(WorkDistrib, recvComputeMap, mapMsg, thisgroup);
}

void WorkDistrib::recvComputeMap(ComputeMapDistribMsg *msg) {
  delete msg;
  DoneMsg *donemsg = new (MsgIndex(DoneMsg)) DoneMsg;
  CSendMsgBranch(WorkDistrib, doneSaveComputeMap, donemsg, thisgroup, 0);
  ComputeMap *computeMap = ComputeMap::Object();

  DebugM(2, "ComputeMap after send!\n");
  for (int i=0; i<computeMap->numComputes(); i++) {
    DebugM(2, "ComputeMap (" << i << ") node = " << computeMap->node(i) << " newNode = " << computeMap->newNode(i) << " type=" << computeMap->type(i) << "\n");
  }
  DebugM(2, "===================================================\n");
}

void WorkDistrib::doneSaveComputeMap(DoneMsg *msg) {
  delete msg;
  DoneMsg *msg2 = new (MsgIndex(DoneMsg)) DoneMsg;
  if (!--saveComputeMapCount) { 
    GeneralSendMsgBranch(saveComputeMapReturnEP, msg2, 0, -1, 
      saveComputeMapReturnChareID);
  }
}


//----------------------------------------------------------------------
// This should only be called on node 0.
//----------------------------------------------------------------------
void WorkDistrib::createHomePatches(void)
{
  int i;
  StringList *current;	//  Pointer used to retrieve configuration items
  Node *node = CLocalBranch(Node,group.node);
  PatchMap *patchMap = PatchMap::Object();
  PatchMgr *patchMgr = CLocalBranch(PatchMgr,group.patchMgr);
  SimParameters *params = node->simParameters;
  Molecule *molecule = node->molecule;
  PDB *pdb = node->pdb;

  int numPatches = patchMap->numPatches();
  int numAtoms = pdb->num_atoms();

  Vector *positions = new Position[numAtoms];
  pdb->get_all_positions(positions);

  Vector *velocities = new Velocity[numAtoms];

  if ( params->initialTemp < 0.0 ) {
    Bool binvels=FALSE;

    //  Reading the veolcities from a PDB
    current = node->configList->find("velocities");

    if (current == NULL) {
      current = node->configList->find("binvelocities");
      binvels = TRUE;
    }

    if (!binvels) {
      velocities_from_PDB(current->data, velocities, numAtoms);
    }
    else {
      velocities_from_binfile(current->data, velocities, numAtoms);
    }
  }
  else {
    // Random velocities for a given temperature
    random_velocities(params->initialTemp, molecule, velocities, numAtoms);
  }

  int numDegFreedom = 3*numAtoms;

  //  If COMMotion == no, remove center of mass motion
  if (!(params->comMove)) {
    remove_com_motion(velocities, molecule, numAtoms);
    numDegFreedom -= 3;
  }

  AtomIDList *atomIDs = new AtomIDList[numPatches];
  PositionList *atomPositions = new PositionList[numPatches];
  VelocityList *atomVelocities = new VelocityList[numPatches];

  Lattice lattice = params->lattice;

  if (params->splitPatch == SPLIT_PATCH_HYDROGEN)
    {
    // split atoms into patched based on helix-group and position
    int aid, pid;
    for(i=0; i < numAtoms; i++)
      {
      if ( ! ( i % 1000 ) )
	{
        DebugM(3,"Assigned " << i << " atoms to patches so far.\n");
        }
      // Assign atoms to patches without splitting hydrogen groups.
      // We know that the hydrogenGroup array is sorted with group parents
      // listed first.  Thus, only change the pid if an atom is a group parent.
      aid = molecule->hydrogenGroup[i].atomID;
      if (molecule->hydrogenGroup[i].isGP)
	pid = patchMap->assignToPatch(positions[aid]);
      // else: don't change pid
      atomIDs[pid].add(aid);
      atomPositions[pid].add(positions[aid]);
      atomVelocities[pid].add(velocities[aid]);
      }
    }
  else
    {
    // split atoms into patched based on position
    for(i=0; i < numAtoms; i++)
      {
      if ( ! ( i % 1000 ) )
	{
	DebugM(3,"Assigned " << i << " atoms to patches so far.\n");
	}
      int pid = patchMap->assignToPatch(positions[i]);
      atomIDs[pid].add(i);
      atomPositions[pid].add(positions[i]);
      atomVelocities[pid].add(velocities[i]);
      }
    }

  delete [] positions;
  delete [] velocities;

  for(i=0; i < numPatches; i++)
  {
    if ( ! ( i % 100 ) )
    {
      DebugM(3,"Created " << i << " patches so far.\n");
    }

    ScaledPosition center(0.5*(patchMap->minX(i)+patchMap->maxX(i)),
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
}

void WorkDistrib::distributeHomePatches() {
  // ref BOC
  Node *node = CLocalBranch(Node,group.node);
  PatchMgr *patchMgr = CLocalBranch(PatchMgr,group.patchMgr);
  // ref singleton
  PatchMap *patchMap = PatchMap::Object();

  // Move patches to the proper node
  for(int i=0;i < patchMap->numPatches(); i++)
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
  delete msg;

  Node *node = CLocalBranch(Node,group.node);

  if (node->myid() != 0)
  {
    DebugM(3,"Saving patch map, compute map\n");
  }
  else
  {
    DebugM(3,"Node 0 patch map built\n");
  }

  mapsArrived = true;
}


//----------------------------------------------------------------------
void WorkDistrib::patchMapInit(void)
{
  PatchMap *patchMap = PatchMap::Object();
  Node *node = CLocalBranch(Node, group.node);
  SimParameters *params = node->simParameters;

  BigReal patchSize = params->patchDimension;

  int xdim, ydim, zdim;
  int xper, yper, zper;
  int xi, yi, zi, pid;
  int i;
  Vector xmin, xmax;
  Vector sysDim, sysMin;

  DebugM(3,"Mapping patches\n");
  node->pdb->find_extremes(&xmin,&xmax);

  DebugM(3,"xmin.x = " << xmin.x << endl);
  DebugM(3,"xmin.y = " << xmin.y << endl);
  DebugM(3,"xmin.z = " << xmin.z << endl);

  DebugM(3,"xmax.x = " << xmax.x << endl);
  DebugM(3,"xmax.y = " << xmax.y << endl);
  DebugM(3,"xmax.z = " << xmax.z << endl);

  Vector center = 0.5 * ( xmax + xmin );

  if ( params->cellBasisVector1.x )
  {
    xper = 1;
    sysDim.x = params->cellBasisVector1.x;
    xdim = (int)(sysDim.x / patchSize);
    if ( xdim < 2 ) xdim = 2;
    iout << iINFO 
	 << "Periodic in x dimension with " << xdim << " patches.\n";
    center.x = params->lattice.origin().x;
  }
  else
  {
    xper = 0;
    sysDim.x = xmax.x - xmin.x;
    xdim = (int)((float)sysDim.x / patchSize);
    if ((xdim * patchSize) < sysDim.x)
      xdim++;
    sysDim.x = xdim * patchSize;
    iout << iINFO
	 << "Non-periodic in x dimension with " << xdim << " patches.\n";
  }

  if ( params->cellBasisVector2.y )
  {
    yper = 1;
    sysDim.y = params->cellBasisVector2.y;
    ydim = (int)((float)sysDim.y / patchSize);
    if ( ydim < 2 ) ydim = 2;
    iout << iINFO
	 << "Periodic in y dimension with " << ydim << " patches.\n";
    center.y = params->lattice.origin().y;
  }
  else
  {
    yper = 0;
    sysDim.y = xmax.y - xmin.y;
    ydim = (int)((float)sysDim.y / patchSize);
    if ((ydim * patchSize) < sysDim.y)
      ydim++;
    sysDim.y = ydim * patchSize;
    iout << iINFO 
	 << "Non-periodic in y dimension with " << ydim << " patches.\n";
  }

  if ( params->cellBasisVector3.z )
  {
    zper = 1;
    sysDim.z = params->cellBasisVector3.z;
    zdim = (int)((float)sysDim.z / patchSize);
    if ( zdim < 2 ) zdim = 2;
    iout << iINFO
	 << "Periodic in z dimension with " << zdim << " patches.\n";
    center.z = params->lattice.origin().z;
  }
  else
  {
    zper = 0;
    sysDim.z = xmax.z - xmin.z;
    zdim = (int)((float)sysDim.z / patchSize);
    if ((zdim * patchSize) < sysDim.z)
      zdim++;
    sysDim.z = zdim * patchSize;
    iout << iINFO
	 << "Non-periodic in z dimension with " << zdim << " patches.\n";
  }

  sysMin = center - 0.5 * sysDim;

  patchMap->setPeriodicity(xper,yper,zper);
  patchMap->allocatePids(xdim, ydim, zdim);

  patchMap->setGridOriginAndLength(sysMin,sysDim);

  Lattice lattice = params->lattice;

  int assignedNode=0;
  for(i=0; i < patchMap->numPatches(); i++)
  {
    pid=patchMap->requestPid(&xi,&yi,&zi); // generates next pid and grid pos
    patchMap->storePatchCoord(pid, 
lattice.scale( Vector(	((float)xi/(float)xdim)*sysDim.x+sysMin.x,
			((float)yi/(float)ydim)*sysDim.y+sysMin.y,
			((float)zi/(float)zdim)*sysDim.z+sysMin.z) ) ,
lattice.scale( Vector(	((float)(xi+1)/(float)xdim)*sysDim.x+sysMin.x,
			((float)(yi+1)/(float)ydim)*sysDim.y+sysMin.y,
			((float)(zi+1)/(float)zdim)*sysDim.z+sysMin.z) ) );
    patchMap->allocateCompute(pid, 100);
    DebugM(3,"Patch " 
        << pid << " is at grid " << xi << " " << yi << " " << zi << ".\n");
  }
}

//----------------------------------------------------------------------
void WorkDistrib::assignNodeToPatch()
{
  int method=1;

  if (method==1)
    assignPatchesRecursiveBisection();
  else
    assignPatchesToLowestLoadNode();

  PatchMap *patchMap = PatchMap::Object();
  int nNodes = Node::Object()->numNodes();
  int *nAtoms = new int[nNodes];
  int numAtoms=0;
  int i;
  for(i=0; i < nNodes; i++)
    nAtoms[i] = 0;

  for(i=0; i < patchMap->numPatches(); i++)
  {
    if (patchMap->patch(i)) {
      numAtoms += patchMap->patch(i)->getNumAtoms();
      nAtoms[patchMap->node(i)] += patchMap->patch(i)->getNumAtoms();
    }
  }

  for(i=0; i < nNodes; i++)
    iout << iINFO 
	 << nAtoms[i] << " atoms assigned to processor " << i << endl;
  iout << iINFO 
       << "Simulation has " << numAtoms << " atoms." << endl;

  delete [] nAtoms;
 
  //  PatchMap::Object()->printPatchMap();
}

//----------------------------------------------------------------------
// void WorkDistrib::assignPatchesSlices() 
// {
//   int pid; 
//   int assignedNode = 0;
//   PatchMap *patchMap = PatchMap::Object();
//   Node *node = CLocalBranch(Node, group.node);

//   int *numAtoms = new int[node->numNodes()];
//   for (int i=0; i<node->numNodes(); i++) {
//     numAtoms[i] = 0;
//   }

//   // Assign patch to node with least atoms assigned.
//   for(pid=0; pid < patchMap->numPatches(); pid++) {
//     assignedNode = 0;
//     for (i=1; i < node->numNodes(); i++) {
//       if (numAtoms[i] < numAtoms[assignedNode]) assignedNode = i;
//     }
//     patchMap->assignNode(pid, assignedNode);
//     numAtoms[assignedNode] += patchMap->patch(pid)->getNumAtoms();

//     /*
//     iout << iINFO << "Patch (" << pid << ") has " 
//       << patchMap->patch(pid)->getNumAtoms() 
//       << " atoms:  Assigned to Node(" << assignedNode << ")\n" 
//       << endi;
//     */
//   }

//   delete[] numAtoms;
// }

//----------------------------------------------------------------------
void WorkDistrib::assignPatchesToLowestLoadNode() 
{
  int pid; 
  int assignedNode = 0;
  PatchMap *patchMap = PatchMap::Object();
  Node *node = CLocalBranch(Node, group.node);

  int *numAtoms = new int[node->numNodes()];
  for (int i=0; i<node->numNodes(); i++) {
    numAtoms[i] = 0;
  }

  // Assign patch to node with least atoms assigned.
  for(pid=0; pid < patchMap->numPatches(); pid++) {
    assignedNode = 0;
    for (i=1; i < node->numNodes(); i++) {
      if (numAtoms[i] < numAtoms[assignedNode]) assignedNode = i;
    }
    patchMap->assignNode(pid, assignedNode);
    numAtoms[assignedNode] += patchMap->patch(pid)->getNumAtoms();

    /*
    iout << iINFO << "Patch (" << pid << ") has " 
      << patchMap->patch(pid)->getNumAtoms() 
      << " atoms:  Assigned to Node(" << assignedNode << ")\n" 
      << endi;
    */
  }

  delete[] numAtoms;
}

//----------------------------------------------------------------------
void WorkDistrib::assignPatchesRecursiveBisection() 
{
  RecBisection recBisec(Node::Object()->numNodes(),PatchMap::Object());
  if ( !recBisec.partition((int *)NULL) ) {
    iout << iWARN 
	 << "WorkDistrib: Recursive bisection fails,"
	 << "invoking least-load algorithm\n";
    assignPatchesToLowestLoadNode();
  }
}

//----------------------------------------------------------------------
void WorkDistrib::mapComputes(void)
{
  PatchMap *patchMap = PatchMap::Object();
  ComputeMap *computeMap = ComputeMap::Object();
  Node *node = CLocalBranch(Node, group.node);

  DebugM(3,"Mapping computes\n");

  // We need to allocate computes for self, 1 and 2 away pairs for
  // electrostatics, and 1 angleForce for each node.  Then I might
  // throw in a few extras, in case I forget some.

  int numPotentialCids = 
    patchMap->numPatches() * (26/2+10) + node->numNodes() * 10;

  computeMap->allocateCids(numPotentialCids);

  // Handle full electrostatics
  if ( node->simParameters->fullDirectOn )
  {
    if ( node->numNodes() > 1 )
      NAMD_die("Full direct electrostatics only works on one processor.");
    else if ( patchMap->xIsPeriodic() ||
		patchMap->yIsPeriodic() || patchMap->zIsPeriodic() )
      NAMD_die("Full direct electrostatics is incompatible with periodic boundary conditions.");
    else
      mapComputeHomePatches(computeFullDirectType);
  }
#ifdef DPMTA
  if ( node->simParameters->FMAOn )
    mapComputeHomePatches(computeDPMTAType);
#endif

  mapComputeNonbonded();
  mapComputeHomePatches(computeNonbondedExclType);
  mapComputeHomePatches(computeBondsType);
  mapComputeHomePatches(computeAnglesType);
  mapComputeHomePatches(computeDihedralsType);
  mapComputeHomePatches(computeImpropersType);

  if ( node->simParameters->sphericalBCOn )
    mapComputePatch(computeSphericalBCType);
  if ( node->simParameters->cylindricalBCOn )
    mapComputePatch(computeCylindricalBCType);
  if ( node->simParameters->constraintsOn )
    mapComputePatch(computeRestraintsType);
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
void WorkDistrib::mapComputePatch(ComputeType type)
{
  PatchMap *patchMap = PatchMap::Object();
  ComputeMap *computeMap = ComputeMap::Object();

  PatchID i;
  ComputeID cid;

  for(i=0; i<patchMap->numPatches(); i++)
  {
    cid=computeMap->storeCompute(patchMap->node(i),1,type);
    computeMap->newPid(cid,i);
    patchMap->newCid(i,cid);
  }

}

//----------------------------------------------------------------------
void WorkDistrib::mapComputeNonbonded(void)
{
  // For each patch, create 1 electrostatic object for self-interaction.
  // Then create 1 for each 1-away and 2-away neighbor which has a larger
  // pid.

  PatchMap *patchMap = PatchMap::Object();
  ComputeMap *computeMap = ComputeMap::Object();
  Node *node = CLocalBranch(Node,group.node);

  PatchID oneAway[PatchMap::MaxOneAway];
  PatchID oneAwayTrans[PatchMap::MaxOneAway];
  PatchID twoAway[PatchMap::MaxTwoAway];

  PatchID i;
  ComputeID cid;
  int numNeighbors;
  int j;
  int numAtoms, numAtoms1, numAtoms2;

  double *pairWork = new double[node->numNodes()];
  for (j=0; j<node->numNodes(); j++) {
    pairWork[j] = 0.;
  }

  const double inspectCost = 0.25;  // Constant cost for checking atoms
  const double computeCost = 1. - inspectCost;
                                    // Cost for actually computing force
  const double weight = inspectCost + computeCost;

  for(i=0; i<patchMap->numPatches(); i++) // do the self 
  {
    // self-interaction
    cid=computeMap->storeCompute(patchMap->node(i),1,computeNonbondedSelfType);
    numAtoms = patchMap->patch(i)->getNumAtoms();
    // pairWork[patchMap->node(i)] += 
      // weight * ((numAtoms*numAtoms)/2. - numAtoms);
    pairWork[patchMap->node(i)] += (numAtoms*numAtoms);
    computeMap->newPid(cid,i);
    patchMap->newCid(i,cid);
  }

  for(i=0; i<patchMap->numPatches(); i++) // do the pairs
  {
    // one-away neighbors
    numNeighbors=patchMap->oneAwayNeighbors(i,oneAway,oneAwayTrans);
    for(j=0;j<numNeighbors;j++)
    {
      if (i < oneAway[j])
      {
        numAtoms1 = patchMap->patch(i)->getNumAtoms();
        numAtoms2 = patchMap->patch(oneAway[j])->getNumAtoms();
	const int distance = 
	  abs(patchMap->xIndex(i)-patchMap->xIndex(oneAway[j])) 
	  + abs(patchMap->yIndex(i)-patchMap->yIndex(oneAway[j])) 
	  + abs(patchMap->zIndex(i)-patchMap->zIndex(oneAway[j]));

        double weight;
        if(distance==1) {
          weight = 0.69;
        } else if(distance==2) {
          weight = 0.32;
        } else if(distance==3) {
          weight = 0.24;
        } else weight = 0;

	if (pairWork[patchMap->node(i)]<pairWork[patchMap->node(oneAway[j])]) {
	    cid=computeMap->storeCompute(patchMap->node(i),2,
				     computeNonbondedPairType);
	    pairWork[patchMap->node(i)] += weight * numAtoms1*numAtoms2;
	} else {
	    cid=computeMap->storeCompute(patchMap->node(oneAway[j]),2,
				     computeNonbondedPairType);
	    pairWork[patchMap->node(oneAway[j])] += weight * numAtoms1*numAtoms2;
	}
	 
	computeMap->newPid(cid,i);
	computeMap->newPid(cid,oneAway[j],oneAwayTrans[j]);
	patchMap->newCid(i,cid);
	patchMap->newCid(oneAway[j],cid);
      }
    }
  }
  /*
  for(i=0; i<node->numNodes(); i++) {
    iout << iINFO << "PairWork on node(" << i << ") = " << pairWork[i] 
      << "\n" << endi;
  }
  */
  delete[] pairWork;
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

//----------------------------------------------------------------------
void WorkDistrib::messageEnqueueWork(Compute *compute) {
  LocalWorkMsg *msg 
    = new (MsgIndex(LocalWorkMsg), sizeof(int)*8) LocalWorkMsg;

  int seq = compute->sequence();
  *CPriorityPtr(msg) = (seq %256) * 256 + compute->priority();

  msg->compute = compute; // pointer is valid since send is to local Pe

  if ( seq < 0 ) {
   CSendMsgBranch(WorkDistrib, enqueueWork, msg, group.workDistrib, CMyPe() );
  }
  else switch ( seq % 3 ) {
  case 0:
   CSendMsgBranch(WorkDistrib, enqueueWorkA, msg, group.workDistrib, CMyPe() );
   break;
  case 1:
   CSendMsgBranch(WorkDistrib, enqueueWorkB, msg, group.workDistrib, CMyPe() );
   break;
  case 2:
   CSendMsgBranch(WorkDistrib, enqueueWorkC, msg, group.workDistrib, CMyPe() );
   break;
  default:
   NAMD_die("WorkDistrib::messageEnqueueWork case statement error!");
   delete msg;
  }
}

void WorkDistrib::enqueueWork(LocalWorkMsg *msg) {
  msg->compute->doWork();
  delete msg;
}

void WorkDistrib::enqueueWorkA(LocalWorkMsg *msg) {
  msg->compute->doWork();
  delete msg;
}

void WorkDistrib::enqueueWorkB(LocalWorkMsg *msg) {
  msg->compute->doWork();
  delete msg;
}

void WorkDistrib::enqueueWorkC(LocalWorkMsg *msg) {
  msg->compute->doWork();
  delete msg;
}

//**********************************************************************
//
//			FUNCTION velocities_from_PDB
//
//   INPUTS:
//      v - Array of vectors to populate
//	filename - name of the PDB filename to read in
//
//	This function reads in a set of initial velocities from a
//      PDB file.  It places the velocities into the array of Vectors
//      passed to it.
//
//***********************************************************************/

void WorkDistrib::velocities_from_PDB(char *filename, 
				      Vector *v, int totalAtoms)
{
  PDB *v_pdb;		//  PDB info from velocity PDB
  int i;

  //  Read the PDB
  v_pdb = new PDB(filename);
  if ( v_pdb == NULL )
  {
    NAMD_die("memory allocation failed in Node::velocities_from_PDB");
  }

  //  Make sure the number of velocities read in matches
  //  the number of atoms we have
  if (v_pdb->num_atoms() != totalAtoms)
  {
    char err_msg[129];

    sprintf(err_msg, "FOUND %d COORDINATES IN VELOCITY PDB!!",
	    v_pdb->num_atoms());

    NAMD_die(err_msg);
  }

  //  Get the entire list of atom info and loop through
  //  them assigning the velocity vector for each one
  v_pdb->get_all_positions(v);

  for (i=0; i<totalAtoms; i++)
  {
    v[i].x *= 0.05;
    v[i].y *= 0.05;
    v[i].z *= 0.05;
  }

  delete v_pdb;
}
//		END OF FUNCTION velocities_from_PDB

//**********************************************************************
//
// 			FUNCTION velocities_from_binfile
//
//    INPUTS:
// 	fname - File name to write velocities to
//	n - Number of atoms in system
//	vels - Array of velocity vectors
//					
//	This function writes out the velocities in binary format.  This is
//     done to preserve accuracy between restarts of namd.
//
//**********************************************************************

void WorkDistrib::velocities_from_binfile(char *fname, Vector *vels, int n)
{
  int32 filen;		//  Number of atoms read from file
  FILE *fp;		//  File descriptor

  //  Open the file and die if the open fails
  if ( (fp = Fopen(fname, "r")) == NULL)
  {
    char errmsg[256];

    sprintf(errmsg, "Unable to open binary velocity file %s", fname);
    NAMD_die(errmsg);
  }

  //  read the number of coordinates in this file
  fread(&filen, sizeof(int32), 1, fp);

  //  Die if this doesn't match the number in our system
  if (filen != n)
  {
    NAMD_die("Number of coordinates in binary velocity file incorrect");
  }

  fread(vels, sizeof(Vector), n, fp);

  Fclose(fp);
}
//               END OF FUNCTION velocities_from_binfile

//**********************************************************************
//
//			FUNCTION random_velocities
//
//   INPUTS:
//	v - array of vectors to populate
//	Temp - Temperature to acheive
//
//	This function assigns a random velocity distribution to a
//   simulation to achieve a desired initial temperature.  The method
//   used here was stolen from the program X-PLOR.
//
//**********************************************************************

void WorkDistrib::random_velocities(BigReal Temp,Molecule *structure,
				    Vector *v, int totalAtoms)
{
  int i, j;		//  Loop counter
  BigReal kbT;		//  Boltzman constant * Temp
  BigReal randnum;	//  Random number from -6.0 to 6.0
  BigReal kbToverM;	//  sqrt(Kb*Temp/Mass)

  kbT = Temp*BOLTZMAN;

  //  Loop through all the atoms and assign velocities in
  //  the x, y and z directions for each one
  for (i=0; i<totalAtoms; i++)
  {
    kbToverM = sqrt(kbT/structure->atommass(i));

    //  The following comment was stolen from X-PLOR where
    //  the following section of code was adapted from.
    
    //  This section generates a Gaussian random
    //  deviate of 0.0 mean and standard deviation RFD for
    //  each of the three spatial dimensions.
    //  The algorithm is a "sum of uniform deviates algorithm"
    //  which may be found in Abramowitz and Stegun,
    //  "Handbook of Mathematical Functions", pg 952.
    for (randnum=0.0, j=0; j<12; j++)
    {
      randnum += NAMD_random();
    }

    randnum -= 6.0;

    v[i].x = randnum*kbToverM;

    for (randnum=0.0, j=0; j<12; j++)
    {
      randnum += NAMD_random();
    }

    randnum -= 6.0;

    v[i].y = randnum*kbToverM;

    for (randnum=0.0, j=0; j<12; j++)
    {
      randnum += NAMD_random();
    }

    randnum -= 6.0;
    
    v[i].z = randnum*kbToverM;
  }
}
/*			END OF FUNCTION random_velocities		*/

//**********************************************************************
//
//			FUNCTION remove_com_motion
//
//   INPUTS:
//	vel - Array of initial velocity vectors
//
//	This function removes the center of mass motion from a molecule.
//
//**********************************************************************

void WorkDistrib::remove_com_motion(Vector *vel, Molecule *structure, int n)
{
  Vector mv;		//  Sum of (mv)_i
  BigReal totalMass=0; 	//  Total mass of system
  int i;			//  Loop counter

  mv.x=0.0;
  mv.y=0.0;
  mv.z=0.0;

  //  Loop through and compute the net momentum
  for (i=0; i<n; i++)
  {
    mv.x += (structure->atommass(i))*vel[i].x;
    mv.y += (structure->atommass(i))*vel[i].y;
    mv.z += (structure->atommass(i))*vel[i].z;
    totalMass += structure->atommass(i);
  }

  mv.x = mv.x/totalMass;
  mv.y = mv.y/totalMass;
  mv.z = mv.z/totalMass;

  //  If any of the velocities really need to change, adjust them
  if ( (fabs(mv.x) > 0.0) || (fabs(mv.y) > 0.0) || (fabs(mv.y) > 0.0) )
  {
    iout << "ADJUSTING COM VELOCITY ("
	 << mv.x << ", " << mv.y << ", " << mv.z  
	 << ") TO REMOVE MOVEMENT\n" << endi;

    for (i=0; i<n; i++)
    {
      vel[i] -= mv;
    }
  }
}
/*			END OF FUNCTION remove_com_motion		*/

#include "WorkDistrib.bot.h"

/***************************************************************************
 * RCS INFORMATION:
 *
 *	$RCSfile: WorkDistrib.C,v $
 *	$Author: brunner $	$Locker:  $		$State: Exp $
 *	$Revision: 1.1035 $	$Date: 1997/09/30 16:22:19 $
 *
 ***************************************************************************
 * REVISION HISTORY:
 *
 * $Log: WorkDistrib.C,v $
 * Revision 1.1035  1997/09/30 16:22:19  brunner
 * Reduced static arrays to allows only 1-away computes to decrease mem usage.
 *
 * Revision 1.1034  1997/09/28 10:19:10  milind
 * Fixed priorities, ReductionMgr etc.
 *
 * Revision 1.1033  1997/08/22 20:12:05  milind
 * Turned on Priorities.
 *
 * Revision 1.1032  1997/08/20 23:27:41  jim
 * Created multiple enqueueWork entry points to aid analysis.
 *
 * Revision 1.1031  97/08/14  15:29:48  15:29:48  brunner (Robert Brunner)
 * More fixes for 32 bit ints in binary restart format
 * 
 * Revision 1.1030  1997/04/22 04:26:03  jim
 * Added atomic restraints (harmonic constraints) via ComputeRestraints class.
 *
 * Revision 1.1029  1997/04/16 22:12:22  brunner
 * Fixed an LdbCoordinator bug, and cleaned up timing and Ldb output some.
 *
 * Revision 1.1028  1997/04/10 15:49:44  brunner
 * Added patch array dimension output and RefineOnly load balancer
 *
 * Revision 1.1027  1997/04/10 14:44:52  milind
 * Changed weights.
 *
 * Revision 1.1026  1997/04/10 09:14:14  ari
 * Final debugging for compute migration / proxy creation for load balancing.
 * Lots of debug code added, mostly turned off now.
 * Fixed bug in PositionBox when Patch had no dependencies.
 * Eliminated use of cout and misuse of iout in numerous places.
 *                                            Ari & Jim
 *
 * Revision 1.1025  1997/04/08 07:09:02  ari
 * Modification for dynamic loadbalancing - moving computes
 * Still bug in new computes or usage of proxies/homepatches.
 * Works if ldbStrategy is none as before.
 *
 * Revision 1.1024  1997/04/07 22:23:31  brunner
 * Changed RB constants so patch distrib equalizes number of atoms, and
 * added weights for initial compute distrib based on self, face neighbor,
 * edge neighbor, or corner neighbor.
 *
 * Revision 1.1023  1997/04/07 21:34:26  brunner
 * Added some information printout about atom distribution.
 *
 * Revision 1.1022  1997/04/07 21:09:58  brunner
 * Added RecBisection for initial patch distrib
 *
 * Revision 1.1021  1997/04/07 14:54:37  nealk
 * Changed fclose() to Fclose() (found in common.[Ch]) to use with popen().
 * Also corrected compilation warnings in Set.[Ch].
 *
 * Revision 1.1020  1997/04/06 22:45:14  ari
 * Add priorities to messages.  Mods to help proxies without computes.
 * Added quick enhancement to end of list insertion of ResizeArray(s)
 *
 * Revision 1.1019  1997/04/03 19:59:13  nealk
 * 1) New Fopen() which handles .Z and .gz files.
 * 2) localWaters and localNonWaters lists on each patch.
 *
 * Revision 1.1018  1997/03/27 17:08:33  nealk
 * Added hydrogen groupings.  Now configuration parameter "splitPatch" determines
 * atom-into-patch distribution.
 *
 * Revision 1.1017  1997/03/27 08:04:27  jim
 * Reworked Lattice to keep center of cell fixed during rescaling.
 *
 * Revision 1.1016  1997/03/19 05:50:14  jim
 * Added ComputeSphericalBC, cleaned up make dependencies.
 *
 * Revision 1.1015  1997/03/15 22:15:35  jim
 * Added ComputeCylindricalBC.  Doesn't break anything but untested and
 * cylinder is along x axis (will fix soon).
 *
 * Revision 1.1014  1997/03/14 21:40:16  ari
 * Reorganized startup to make possible inital load
 * balancing by changing methods in WorkDistrib.
 * Also made startup more transparent and easier
 * to modify.
 *
 * Revision 1.1012  1997/03/10 17:40:18  ari
 * UniqueSet changes - some more commenting and cleanup
 *
 * Revision 1.1011  1997/03/08 21:09:04  jim
 * Patches are now centered on system so FMA results should match NAMD 1.X.
 *
 * Revision 1.1010  1997/03/06 22:18:16  brunner
 * Made utility functions private functions of WorkDistrib.
 *
 * Revision 1.1009  1997/03/06 22:06:10  ari
 * Removed Compute.ci
 * Comments added - more code cleaning
 *
 * Revision 1.1008  1997/02/28 16:13:55  nealk
 * Turned off debugging code.
 *
 * Revision 1.1007  1997/02/26 23:18:45  jim
 * Now should read binary coordinate files - untested.
 *
 * Revision 1.1006  1997/02/26 16:53:19  ari
 * Cleaning and debuging for memory leaks.
 * Adding comments.
 * Removed some dead code due to use of Quiescense detection.
 *
 * Revision 1.1005  1997/02/14 20:24:57  jim
 * Patches are now sized according to config file.
 *
 * Revision 1.1004  1997/02/13 22:27:04  jim
 * Added inital velocity code from NAMD 1.
 * Reading velocity pdb file appears to work.
 * Reading binary velociy file should work but is untested.
 * Random velocites appears to work but differs from NAMD 1.
 *
 * Revision 1.1003  1997/02/13 04:43:17  jim
 * Fixed initial hanging (bug in PatchMap, but it still shouldn't have
 * happened) and saved migration messages in the buffer from being
 * deleted, but migration still dies (even on one node).
 *
 * Revision 1.1002  1997/02/11 18:51:58  ari
 * Modified with #ifdef DPMTA to safely eliminate DPMTA codes
 * fixed non-buffering of migration msgs
 * Migration works on multiple processors
 *
 * Revision 1.1001  1997/02/07 17:39:42  ari
 * More debugging for atomMigration.
 * Using -w on CC got us some minor fixes
 * using purify got us a major memory problem due to bad sizing of dummy force
 *
 * Revision 1.1000  1997/02/06 15:59:26  ari
 * Resetting CVS to merge branches back into the main trunk.
 * We will stick to main trunk development as suggested by CVS manual.
 * We will set up tags to track fixed points of development/release
 * as suggested by CVS manual - all praise the CVS manual.
 *
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

