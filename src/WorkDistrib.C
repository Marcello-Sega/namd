/**
***  Copyright (c) 1995, 1996, 1997, 1998, 1999, 2000 by
***  The Board of Trustees of the University of Illinois.
***  All rights reserved.
**/

/*
   Currently, WorkDistrib generates the layout of the Patches,
   directs the construction and distribution of Computes and
   associates Computes with Patches.
*/

#include <stdio.h>

#include "charm++.h"

#include "ProcessorPrivate.h"

#include "BOCgroup.h"
#include "WorkDistrib.decl.h"
#include "WorkDistrib.h"

#include "main.decl.h"
#include "main.h"
#include "Node.h"
#include "PatchMgr.h"
#include "PatchMap.inl"
#include "NamdTypes.h"
#include "PDB.h"
#include "SimParameters.h"
#include "Molecule.h"
#include "NamdOneTools.h"
#include "Compute.h"
#include "ComputeMap.h"
#include "RecBisection.h"
#include "Random.h"

//#define MIN_DEBUG_LEVEL 4
//#define DEBUGM
#include "Debug.h"

//======================================================================
// Public functions
//----------------------------------------------------------------------
WorkDistrib::WorkDistrib()
{
  CpvAccess(BOCclass_group).workDistrib = thisgroup;
  mapsArrived = false;
  awaitingMaps = false;
}

//----------------------------------------------------------------------
WorkDistrib::~WorkDistrib(void)
{ }


//----------------------------------------------------------------------
void WorkDistrib::sendMaps(void)
{
  MapDistribMsg *mapMsg = new MapDistribMsg ;

  mapMsg->patchMap = PatchMap::Object();
  mapMsg->computeMap = ComputeMap::Object();

  CProxy_WorkDistrib(thisgroup).saveMaps(mapMsg);
  mapsArrived = true;
}

//----------------------------------------------------------------------
void WorkDistrib::saveComputeMapChanges(int ep, int chareID)
{
  saveComputeMapReturnEP = ep;
  saveComputeMapReturnChareID = chareID;
  saveComputeMapCount = CkNumPes();

  ComputeMap *computeMap = ComputeMap::Object();

  int i;
  for (i=0; i<computeMap->numComputes(); i++) {
    DebugM(3, "ComputeMap (" << i << ") node = " << computeMap->node(i) << " newNode = " << computeMap->newNode(i) << "\n");
  }
  
  ComputeMapChangeMsg *mapMsg 
    = new ComputeMapChangeMsg ;

  mapMsg->numNewNodes = computeMap->numComputes();
  for(i=0; i<computeMap->numComputes(); i++)
    mapMsg->newNodes[i] = computeMap->newNode(i);

  CProxy_WorkDistrib(thisgroup).recvComputeMapChanges(mapMsg);
}

void WorkDistrib::recvComputeMapChanges(ComputeMapChangeMsg *msg) {
  
  ComputeMap *computeMap = ComputeMap::Object();
  int i;
  for(i=0; i<computeMap->numComputes(); i++)
    computeMap->setNewNode(i,msg->newNodes[i]);

  delete msg;

  CProxy_WorkDistrib(thisgroup).doneSaveComputeMap(0);

  DebugM(2, "ComputeMap after send!\n");
  for (i=0; i<computeMap->numComputes(); i++) {
    DebugM(2, "ComputeMap (" << i << ") node = " << computeMap->node(i) << " newNode = " << computeMap->newNode(i) << " type=" << computeMap->type(i) << "\n");
  }
  DebugM(2, "===================================================\n");
}

void WorkDistrib::doneSaveComputeMap() {
  if (!--saveComputeMapCount) { 
    CkSendMsgBranch(saveComputeMapReturnEP, CkAllocMsg(0,0,0), 0, saveComputeMapReturnChareID);
  }
}


//----------------------------------------------------------------------
// This should only be called on node 0.
//----------------------------------------------------------------------
void WorkDistrib::createHomePatches(void)
{
  int i;
  StringList *current;	//  Pointer used to retrieve configuration items
  CProxy_Node nd(CpvAccess(BOCclass_group).node);
  Node *node = nd.ckLocalBranch();
  PatchMap *patchMap = PatchMap::Object();
  CProxy_PatchMgr pm(CpvAccess(BOCclass_group).patchMgr);
  PatchMgr *patchMgr = pm.ckLocalBranch();
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

  FullAtomList *atoms = new FullAtomList[numPatches];

  const Lattice lattice = params->lattice;

  if (params->splitPatch == SPLIT_PATCH_HYDROGEN)
    {
    // split atoms into patched based on helix-group and position
    int aid, pid=0;
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
	pid = patchMap->assignToPatch(positions[aid],lattice);
      // else: don't change pid
      FullAtom a;
      a.id = aid;
      a.position = positions[aid];
      a.velocity = velocities[aid];
      atoms[pid].add(a);
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
      int pid = patchMap->assignToPatch(positions[i],lattice);
      FullAtom a;
      a.id = i;
      a.position = positions[i];
      a.velocity = velocities[i];
      atoms[pid].add(a);
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

    ScaledPosition center(0.5*(patchMap->min_a(i)+patchMap->max_a(i)),
			  0.5*(patchMap->min_b(i)+patchMap->max_b(i)),
			  0.5*(patchMap->min_c(i)+patchMap->max_c(i)));

    int n = atoms[i].size();
    FullAtom *a = atoms[i].begin();
    int j;
    for(j=0; j < n; j++)
    {
      int aid = a[j].id;
      a[j].position = lattice.nearest(
		a[j].position, center, &(a[j].transform));
      a[j].mass = molecule->atommass(aid);
      a[j].charge = molecule->atomcharge(aid);
      if (params->splitPatch == SPLIT_PATCH_HYDROGEN) {
        if ( molecule->is_hydrogenGroupParent(aid) ) {
          a[j].hydrogenGroupSize = molecule->get_groupSize(aid);
        } else {
          a[j].hydrogenGroupSize = 0;
        }
      } else {
        a[j].hydrogenGroupSize = 1;
      }
      a[j].nonbondedGroupSize = a[j].hydrogenGroupSize;
      a[j].atomFixed = molecule->is_atom_fixed(aid) ? 1 : 0;
    }

    int size, allfixed, k;
    for(j=0; j < n; j+=size) {
      size = a[j].hydrogenGroupSize;
      if ( ! size ) {
        NAMD_bug("Mother atom with hydrogenGroupSize of 0!");
      }
      allfixed = 1;
      for ( k = 0; k < size; ++k ) {
        allfixed = ( allfixed && (a[j+k].atomFixed) );
      }
      for ( k = 0; k < size; ++k ) {
        a[j+k].groupFixed = allfixed ? 1 : 0;
      }
    }

    patchMgr->createHomePatch(i,atoms[i]);
  }

  delete [] atoms;
}

void WorkDistrib::distributeHomePatches() {
  // ref BOC
  CProxy_Node nd(CpvAccess(BOCclass_group).node);
  Node *node = nd.ckLocalBranch();
  CProxy_PatchMgr pm(CpvAccess(BOCclass_group).patchMgr);
  PatchMgr *patchMgr = pm.ckLocalBranch();
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

  CProxy_Node nd(CpvAccess(BOCclass_group).node);
  Node *node = nd.ckLocalBranch();

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
  CProxy_Node nd(CpvAccess(BOCclass_group).node);
  Node *node = nd.ckLocalBranch();
  SimParameters *params = node->simParameters;
  Lattice lattice = params->lattice;

  BigReal patchSize = params->patchDimension;

  int xdim, ydim, zdim;
  int xper, yper, zper;
  int xi, yi, zi, pid;
  int i;
  ScaledPosition xmin, xmax;
  ScaledPosition sysDim, sysMin;

  DebugM(3,"Mapping patches\n");
  // Need to use full box for FMA to match NAMD 1.X results.
  if ( params->FMAOn ) {
    node->pdb->find_extremes(&(xmin.x),&(xmax.x),lattice.a_r());
    node->pdb->find_extremes(&(xmin.y),&(xmax.y),lattice.b_r());
    node->pdb->find_extremes(&(xmin.z),&(xmax.z),lattice.c_r());
  // Otherwise, this allows a small number of stray atoms.
  } else {
    node->pdb->find_extremes(&(xmin.x),&(xmax.x),lattice.a_r(),0.99);
    node->pdb->find_extremes(&(xmin.y),&(xmax.y),lattice.b_r(),0.99);
    node->pdb->find_extremes(&(xmin.z),&(xmax.z),lattice.c_r(),0.99);
  }

  patchMap->initialize(xmin,xmax,lattice,patchSize);

/*
  ScaledPosition center = 0.5 * ( xmax + xmin );

  iout << iINFO << "PATCH GRID IS ";

  if ( lattice.a_p() )
  {
    xper = 1;
    sysDim.x = cross(lattice.b(),lattice.c()).unit() * lattice.a();
    xdim = (int)(sysDim.x / patchSize);
    if ( xdim < 2 ) xdim = 2;
    iout << xdim << " (PERIODIC)";
    center.x = 0.0;
  }
  else
  {
    xper = 0;
    sysDim.x = xmax.x - xmin.x;
    xdim = (int)((float)sysDim.x / patchSize);
    if ((xdim * patchSize) < sysDim.x)
      xdim++;
    sysDim.x = xdim * patchSize;
    iout << xdim;
  }

  iout << " BY ";

  if ( lattice.b_p() )
  {
    yper = 1;
    sysDim.y = cross(lattice.c(),lattice.a()).unit() * lattice.b();
    ydim = (int)((float)sysDim.y / patchSize);
    if ( ydim < 2 ) ydim = 2;
    iout << ydim << " (PERIODIC)";
    center.y = 0.0;
  }
  else
  {
    yper = 0;
    sysDim.y = xmax.y - xmin.y;
    ydim = (int)((float)sysDim.y / patchSize);
    if ((ydim * patchSize) < sysDim.y)
      ydim++;
    sysDim.y = ydim * patchSize;
    iout << ydim;
  }

  iout << " BY ";

  if ( lattice.c_p() )
  {
    zper = 1;
    sysDim.z = cross(lattice.a(),lattice.b()).unit() * lattice.c();
    zdim = (int)((float)sysDim.z / patchSize);
    if ( zdim < 2 ) zdim = 2;
    iout << zdim << " (PERIODIC)";
    center.z = 0.0;
  }
  else
  {
    zper = 0;
    sysDim.z = xmax.z - xmin.z;
    zdim = (int)((float)sysDim.z / patchSize);
    if ((zdim * patchSize) < sysDim.z)
      zdim++;
    sysDim.z = zdim * patchSize;
    iout << zdim;
  }

  iout << "\n" << endi;

  sysMin = center - 0.5 * sysDim;

  patchMap->setPeriodicity(xper,yper,zper);
  patchMap->allocatePids(xdim, ydim, zdim);

  patchMap->setGridOriginAndLength(sysMin,sysDim);

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
*/
}

//----------------------------------------------------------------------
void WorkDistrib::assignNodeToPatch()
{
  int method=1;

  PatchMap *patchMap = PatchMap::Object();
  int nNodes = Node::Object()->numNodes();
  if (nNodes > patchMap->numPatches())
    //NAMD_die("Man, Tiny Elvis, WorkDistrib::assignNodeToPatch not enough patches!");
    assignPatchesRoundRobin();
  else if (nNodes == patchMap->numPatches())
    assignPatchesRoundRobin();
  else if (method==1)
    assignPatchesRecursiveBisection();
  else
    assignPatchesToLowestLoadNode();

  int *nAtoms = new int[nNodes];
  int numAtoms=0;
  int i;
  for(i=0; i < nNodes; i++)
    nAtoms[i] = 0;

  for(i=0; i < patchMap->numPatches(); i++)
  {
    //    iout << iINFO << "Patch " << i << " has " 
    //	 << patchMap->patch(i)->getNumAtoms() << " atoms and "
    //	 << patchMap->patch(i)->getNumAtoms() * 
    //            patchMap->patch(i)->getNumAtoms() 
    //	 << " pairs.\n" << endi;

    if (patchMap->patch(i)) {
      numAtoms += patchMap->patch(i)->getNumAtoms();
      nAtoms[patchMap->node(i)] += patchMap->patch(i)->getNumAtoms();
    }
  }

//  for(i=0; i < nNodes; i++)
//    iout << iINFO 
//	 << nAtoms[i] << " atoms assigned to node " << i << "\n" << endi;
  if ( numAtoms != Node::Object()->molecule->numAtoms ) {
    NAMD_die("Incorrect atom count in WorkDistrib::assignNodeToPatch\n");
  }

  delete [] nAtoms;
 
  //  PatchMap::Object()->printPatchMap();
}

//----------------------------------------------------------------------
// void WorkDistrib::assignPatchesSlices() 
// {
//   int pid; 
//   int assignedNode = 0;
//   PatchMap *patchMap = PatchMap::Object();
//   Node *node = CLocalBranch(Node, CpvAccess(BOCclass_group).node);

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
  CProxy_Node nd(CpvAccess(BOCclass_group).node);
  Node *node = nd.ckLocalBranch();

  int *load = new int[node->numNodes()];
  for (int i=0; i<node->numNodes(); i++) {
    load[i] = 0;
  }

  // Assign patch to node with least atoms assigned.
  for(pid=0; pid < patchMap->numPatches(); pid++) {
    assignedNode = 0;
    for (int i=1; i < node->numNodes(); i++) {
      if (load[i] < load[assignedNode]) assignedNode = i;
    }
    patchMap->assignNode(pid, assignedNode);
    load[assignedNode] += patchMap->patch(pid)->getNumAtoms() + 1;

    /*
    iout << iINFO << "Patch (" << pid << ") has " 
      << patchMap->patch(pid)->getNumAtoms() 
      << " atoms:  Assigned to Node(" << assignedNode << ")\n" 
      << endi;
    */
  }

  delete[] load;
}

//----------------------------------------------------------------------
void WorkDistrib::assignPatchesRoundRobin() 
{
  int pid; 
  int assignedNode = 0;
  PatchMap *patchMap = PatchMap::Object();
  CProxy_Node nd(CpvAccess(BOCclass_group).node);
  Node *node = nd.ckLocalBranch();

  for(pid=0; pid < patchMap->numPatches(); pid++) {
    assignedNode = pid % node->numNodes();
    patchMap->assignNode(pid, assignedNode);
  }
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
  CProxy_Node nd(CpvAccess(BOCclass_group).node);
  Node *node = nd.ckLocalBranch();

  DebugM(3,"Mapping computes\n");

  // We need to allocate computes for self, 1 and 2 away pairs for
  // electrostatics, and 1 angleForce for each node.  Then I might
  // throw in a few extras, in case I forget some.

#define MAX_SELF_PARTITIONS 100
#define MAX_PAIR_PARTITIONS 10

/*
  int numPotentialCids =
	patchMap->numPatches() *
		(13 * MAX_PAIR_PARTITIONS + MAX_SELF_PARTITIONS + 10) +
	node->numNodes() * 20;
*/
  int numPotentialCids =
	patchMap->numPatches() *
		(13 * node->simParameters->maxPairPart + 
		node->simParameters->maxSelfPart + 10) +
	node->numNodes() * 20;
  //iout << iINFO << "numPotentialCids: " << numPotentialCids << "\n" << endl;
  computeMap->allocateCids(numPotentialCids);

  // Handle full electrostatics
  if ( node->simParameters->fullDirectOn )
    mapComputeHomePatches(computeFullDirectType);
  if ( node->simParameters->FMAOn )
#ifdef DPMTA
    mapComputeHomePatches(computeDPMTAType);
#else
    NAMD_die("This binary does not include DPMTA (FMA).");
#endif
  if ( node->simParameters->PMEOn ) {
#ifdef DPME
    if ( node->simParameters->useDPME )
      mapComputeHomePatches(computeDPMEType);
    else {
      mapComputeHomePatches(computePmeType);
    }
#else
    mapComputeHomePatches(computePmeType);
#endif
  }

  if ( node->simParameters->globalForcesOn )
    mapComputeHomePatches(computeGlobalType);

  mapComputeNonbonded();
  mapComputeHomePatches(computeBondsType);
  mapComputeHomePatches(computeAnglesType);
  mapComputeHomePatches(computeDihedralsType);
  mapComputeHomePatches(computeImpropersType);
  mapComputePatch(computeSelfBondsType);
  mapComputePatch(computeSelfAnglesType);
  mapComputePatch(computeSelfDihedralsType);
  mapComputePatch(computeSelfImpropersType);

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
  CProxy_Node nd(CpvAccess(BOCclass_group).node);
  Node *node = nd.ckLocalBranch();

  int i;

  int numNodes = node->numNodes();

  if (numNodes > patchMap->numPatches()) numNodes=patchMap->numPatches();
  ComputeID *cid = new ComputeID[numNodes];

  for(i=0; i<numNodes; i++)
  {
    cid[i]=computeMap->storeCompute(i,patchMap->numPatches(),type);
  }

  PatchID j;

  for(j=0;j<patchMap->numPatches();j++)
  {
    patchMap->newCid(j,cid[patchMap->node(j)]);
    computeMap->newPid(cid[patchMap->node(j)],j);
  }

  delete [] cid;
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

#ifndef WIN32
static inline int max(int x1,int x2)
{
  return (x1 > x2) ? x1 : x2;
}
#endif

//----------------------------------------------------------------------
void WorkDistrib::mapComputeNonbonded(void)
{
  // For each patch, create 1 electrostatic object for self-interaction.
  // Then create 1 for each 1-away and 2-away neighbor which has a larger
  // pid.

  PatchMap *patchMap = PatchMap::Object();
  ComputeMap *computeMap = ComputeMap::Object();
  CProxy_Node nd(CpvAccess(BOCclass_group).node);
  Node *node = nd.ckLocalBranch();

  PatchID oneAway[PatchMap::MaxOneAway];
  PatchID oneAwayTrans[PatchMap::MaxOneAway];
  // PatchID twoAway[PatchMap::MaxTwoAway];

  PatchID i;
  ComputeID cid;
  int numNeighbors;
  int j;

  double *pairWork = new double[node->numNodes()];
  for (j=0; j<node->numNodes(); j++) {
    pairWork[j] = 0.;
  }

  // const double inspectCost = 0.25;  // Constant cost for checking atoms
  // const double computeCost = 1. - inspectCost;
                                    // Cost for actually computing force
  // const double weight = inspectCost + computeCost;

  for(i=0; i<patchMap->numPatches(); i++) // do the self 
  {
    int numAtoms = patchMap->patch(i)->getNumAtoms();
    int numPartitions;
    if (node->simParameters->numAtomsSelf == 0) {
      numPartitions = 1 + (numAtoms > 50) + (numAtoms*numAtoms)/50000;
    }
    else {
      numPartitions =
        numAtoms*numAtoms / (double)(node->simParameters->numAtomsSelf*node->simParameters->numAtomsSelf) + 0.5;
      if (numPartitions < 1) numPartitions = 1;
    }
    if ( numPartitions > node->simParameters->maxSelfPart )
			numPartitions = node->simParameters->maxSelfPart;
    // self-interaction
    DebugM(4,"Mapping " << numPartitions << " ComputeNonbondedSelf objects for patch " << i << "\n");
//    iout <<"Self numPartitions = " <<numPartitions <<" numAtoms " <<numAtoms <<endl;
    for(int partition=0; partition < numPartitions; partition++)
    {
      cid=computeMap->storeCompute(patchMap->node(i),1,
				   computeNonbondedSelfType,
				   partition,numPartitions);
      // pairWork[patchMap->node(i)] += 
      // weight * ((numAtoms*numAtoms)/2. - numAtoms);
      pairWork[patchMap->node(i)] += (numAtoms*numAtoms) / numPartitions;
      computeMap->newPid(cid,i);
      patchMap->newCid(i,cid);
    }
  }

//   for(i=0; i<patchMap->numPatches(); i++) // do the pairs
//   {
//     // one-away neighbors
//     numNeighbors=patchMap->oneAwayNeighbors(i,oneAway,oneAwayTrans);
//     for(j=0;j<numNeighbors;j++)
//     {
//       if (i < oneAway[j])
//       {
//         numAtoms1 = patchMap->patch(i)->getNumAtoms();
//         numAtoms2 = patchMap->patch(oneAway[j])->getNumAtoms();
// 	const int distance = 
// 	  abs(patchMap->index_a(i)-patchMap->index_a(oneAway[j])) 
// 	  + abs(patchMap->index_b(i)-patchMap->index_b(oneAway[j])) 
// 	  + abs(patchMap->index_c(i)-patchMap->index_c(oneAway[j]));

//         double weight;
//         if(distance==1) {
//           weight = 0.69;
//         } else if(distance==2) {
//           weight = 0.32;
//         } else if(distance==3) {
//           weight = 0.24;
//         } else weight = 0;

// 	if (pairWork[patchMap->node(i)]<pairWork[patchMap->node(oneAway[j])]) {
// 	    cid=computeMap->storeCompute(patchMap->node(i),2,
// 				     computeNonbondedPairType);
// 	    pairWork[patchMap->node(i)] += weight * numAtoms1*numAtoms2;
// 	} else {
// 	    cid=computeMap->storeCompute(patchMap->node(oneAway[j]),2,
// 				     computeNonbondedPairType);
// 	    pairWork[patchMap->node(oneAway[j])] += weight * numAtoms1*numAtoms2;
// 	}
	 
// 	computeMap->newPid(cid,i);
// 	computeMap->newPid(cid,oneAway[j],oneAwayTrans[j]);
// 	patchMap->newCid(i,cid);
// 	patchMap->newCid(oneAway[j],cid);
//       }
//     }
//   }

  for(int p1=0; p1 <patchMap->numPatches(); p1++) // do the pairs
  {
    // one-away neighbors
    numNeighbors=patchMap->oneAwayNeighbors(p1,oneAway,oneAwayTrans);
    for(j=0;j<numNeighbors;j++)
    {
      if (p1 < oneAway[j])
      {
	int p2 = oneAway[j];
	int numAtoms1 = patchMap->patch(p1)->getNumAtoms();
	int numAtoms2 = patchMap->patch(p2)->getNumAtoms();
	const int distance =
 	  ( patchMap->index_a(p1) == patchMap->index_a(p2) ? 0 : 1 ) +
 	  ( patchMap->index_b(p1) == patchMap->index_b(p2) ? 0 : 1 ) +
 	  ( patchMap->index_c(p1) == patchMap->index_c(p2) ? 0 : 1 );
        int numPartitions;
	int divide = 0;
        if (distance <= 1) {
	  divide = node->simParameters->numAtomsPair;
	}
	else if (distance == 2) {
	  divide = node->simParameters->numAtomsPair2;
	}
	else {
	  divide = 0;
	}
	if (divide == 0) {
          numPartitions = 1 + (numAtoms1*numAtoms2 > 2500) + (numAtoms1*numAtoms2)/100000;
	}
	else {
          numPartitions = numAtoms1*numAtoms2/(double)(divide*divide) + 0.5;
          if ( numPartitions < 1 ) numPartitions = 1;
	}
        if ( numPartitions > node->simParameters->maxPairPart )
			numPartitions = node->simParameters->maxPairPart;
//	if ( numPartitions > 1 ) iout << "Mapping " << numPartitions << " ComputeNonbondedPair objects for patches " << p1 << "(" << numAtoms1 << ") and " << p2 << "(" << numAtoms2 << ")\n" << endi;
	for(int partition=0; partition < numPartitions; partition++)
	{
	  cid=computeMap->storeCompute(
		patchMap->node(patchMap->downstream(p1,p2)),
		2,computeNonbondedPairType,partition,numPartitions);
	  computeMap->newPid(cid,p1);
	  computeMap->newPid(cid,p2,oneAwayTrans[j]);
	  patchMap->newCid(p1,cid);
	  patchMap->newCid(p2,cid);
        }

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
  LocalWorkMsg *msg = compute->localWorkMsg;
  CkSetQueueing(msg, CK_QUEUEING_IFIFO);
  int seq = compute->sequence();

  if ( seq < 0 ) {
    *((int*) CkPriorityPtr(msg)) = 128 + compute->priority();
  } else {
    *((int*) CkPriorityPtr(msg)) = 128 + (seq %256) * 256 + compute->priority();
  }

  msg->compute = compute; // pointer is valid since send is to local Pe
  int type = compute->type();

  CProxy_WorkDistrib wdProxy(CpvAccess(BOCclass_group).workDistrib);
  switch ( type ) {
  case computeBondsType:
  case computeSelfBondsType:
    wdProxy.enqueueBonds(msg,CkMyPe());
    break;
  case computeAnglesType:
  case computeSelfAnglesType:
    wdProxy.enqueueAngles(msg,CkMyPe());
    break;
  case computeDihedralsType:
  case computeSelfDihedralsType:
    wdProxy.enqueueDihedrals(msg,CkMyPe());
    break;
  case computeImpropersType:
  case computeSelfImpropersType:
    wdProxy.enqueueImpropers(msg,CkMyPe());
    break;
  case computeNonbondedSelfType:
    switch ( seq % 2 ) {
    case 0:
      wdProxy.enqueueSelfA(msg,CkMyPe());
      break;
    case 1:
      wdProxy.enqueueSelfB(msg,CkMyPe());
      break;
    default:
      NAMD_bug("WorkDistrib::messageEnqueueSelf case statement error!");
    }
    break;
  case computeNonbondedPairType:
    switch ( seq % 2 ) {
    case 0:
      wdProxy.enqueueWorkA(msg,CkMyPe());
      break;
    case 1:
      wdProxy.enqueueWorkB(msg,CkMyPe());
      break;
    case 2:
      wdProxy.enqueueWorkC(msg,CkMyPe());
      break;
    default:
      NAMD_bug("WorkDistrib::messageEnqueueWork case statement error!");
    }
    break;
  case computePmeType:
    wdProxy.enqueuePme(msg,CkMyPe());
    break;
  default:
    wdProxy.enqueueWork(msg,CkMyPe());
  }
}

void WorkDistrib::enqueueWork(LocalWorkMsg *msg) {
  msg->compute->doWork();
  if ( msg->compute->localWorkMsg != msg )
    NAMD_bug("WorkDistrib LocalWorkMsg recycling failed!");
}

void WorkDistrib::enqueueBonds(LocalWorkMsg *msg) {
  msg->compute->doWork();
  if ( msg->compute->localWorkMsg != msg )
    NAMD_bug("WorkDistrib LocalWorkMsg recycling failed!");
}

void WorkDistrib::enqueueAngles(LocalWorkMsg *msg) {
  msg->compute->doWork();
  if ( msg->compute->localWorkMsg != msg )
    NAMD_bug("WorkDistrib LocalWorkMsg recycling failed!");
}

void WorkDistrib::enqueueDihedrals(LocalWorkMsg *msg) {
  msg->compute->doWork();
  if ( msg->compute->localWorkMsg != msg )
    NAMD_bug("WorkDistrib LocalWorkMsg recycling failed!");
}

void WorkDistrib::enqueueImpropers(LocalWorkMsg *msg) {
  msg->compute->doWork();
  if ( msg->compute->localWorkMsg != msg )
    NAMD_bug("WorkDistrib LocalWorkMsg recycling failed!");
}

void WorkDistrib::enqueuePme(LocalWorkMsg *msg) {
  msg->compute->doWork();
  if ( msg->compute->localWorkMsg != msg )
    NAMD_bug("WorkDistrib LocalWorkMsg recycling failed!");
}

void WorkDistrib::enqueueSelfA(LocalWorkMsg *msg) {
  msg->compute->doWork();
  if ( msg->compute->localWorkMsg != msg )
    NAMD_bug("WorkDistrib LocalWorkMsg recycling failed!");
}

void WorkDistrib::enqueueSelfB(LocalWorkMsg *msg) {
  msg->compute->doWork();
  if ( msg->compute->localWorkMsg != msg )
    NAMD_bug("WorkDistrib LocalWorkMsg recycling failed!");
}

void WorkDistrib::enqueueWorkA(LocalWorkMsg *msg) {
  msg->compute->doWork();
  if ( msg->compute->localWorkMsg != msg )
    NAMD_bug("WorkDistrib LocalWorkMsg recycling failed!");
}

void WorkDistrib::enqueueWorkB(LocalWorkMsg *msg) {
  msg->compute->doWork();
  if ( msg->compute->localWorkMsg != msg )
    NAMD_bug("WorkDistrib LocalWorkMsg recycling failed!");
}

void WorkDistrib::enqueueWorkC(LocalWorkMsg *msg) {
  msg->compute->doWork();
  if ( msg->compute->localWorkMsg != msg )
    NAMD_bug("WorkDistrib LocalWorkMsg recycling failed!");
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
    v[i].x *= PDBVELINVFACTOR;
    v[i].y *= PDBVELINVFACTOR;
    v[i].z *= PDBVELINVFACTOR;
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
  read_binary_file(fname,vels,n);
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
  Random vel_random(Node::Object()->simParameters->randomSeed);

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
      randnum += vel_random.uniform();
    }

    randnum -= 6.0;

    v[i].x = randnum*kbToverM;

    for (randnum=0.0, j=0; j<12; j++)
    {
      randnum += vel_random.uniform();
    }

    randnum -= 6.0;

    v[i].y = randnum*kbToverM;

    for (randnum=0.0, j=0; j<12; j++)
    {
      randnum += vel_random.uniform();
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
  Vector mv(0,0,0);		//  Sum of (mv)_i
  BigReal totalMass=0; 	//  Total mass of system
  int i;			//  Loop counter

  //  Loop through and compute the net momentum
  for (i=0; i<n; i++)
  {
    BigReal mass = structure->atommass(i);
    mv += mass * vel[i];
    totalMass += mass;
  }

  mv /= totalMass;

  iout << iINFO << "REMOVING COM VELOCITY "
	<< ( PDBVELFACTOR * mv ) << "\n" << endi;

  for (i=0; i<n; i++) { vel[i] -= mv; }

}
/*			END OF FUNCTION remove_com_motion		*/

#include "WorkDistrib.def.h"

