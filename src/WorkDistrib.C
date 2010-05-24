/**
***  Copyright (c) 1995, 1996, 1997, 1998, 1999, 2000 by
***  The Board of Trustees of the University of Illinois.
***  All rights reserved.
**/

/*****************************************************************************
 * $Source: /home/cvs/namd/cvsroot/namd2/src/WorkDistrib.C,v $
 * $Author: gzheng $
 * $Date: 2010/05/24 13:52:37 $
 * $Revision: 1.1206 $
 *****************************************************************************/

/** \file WorkDistrib.C
 *  Currently, WorkDistrib generates the layout of the Patches,
 *  directs the construction and distribution of Computes and
 *  associates Computes with Patches.
 */

#include <stdio.h>

#include "InfoStream.h"
#include "ProcessorPrivate.h"
#include "BOCgroup.h"
#include "WorkDistrib.decl.h"
#include "WorkDistrib.h"
#include "Lattice.h"
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
#include "varsizemsg.h"
#include "ProxyMgr.h"
#include "Priorities.h"

//#define DEBUGM
#define MIN_DEBUG_LEVEL 2
#include "Debug.h"


class ComputeMapChangeMsg : public CMessage_ComputeMapChangeMsg
{
public:

  int numNewNodes;
  int *newNodes;

//  VARSIZE_DECL(ComputeMapChangeMsg);
};

/*
VARSIZE_MSG(ComputeMapChangeMsg,
  VARSIZE_ARRAY(newNodes);
)
*/


//======================================================================
// Public functions
//----------------------------------------------------------------------
WorkDistrib::WorkDistrib()
{
  CkpvAccess(BOCclass_group).workDistrib = thisgroup;
  mapsArrived = false;
  awaitingMaps = false;
}

//----------------------------------------------------------------------
WorkDistrib::~WorkDistrib(void)
{ }


//----------------------------------------------------------------------
void WorkDistrib::saveComputeMapChanges(int ep, CkGroupID chareID)
{
  saveComputeMapReturnEP = ep;
  saveComputeMapReturnChareID = chareID;
  saveComputeMapCount = CkNumPes();

  ComputeMap *computeMap = ComputeMap::Object();

  int i;
  int nc = computeMap->numComputes();
  
  ComputeMapChangeMsg *mapMsg = new (nc, 0) ComputeMapChangeMsg ;

  mapMsg->numNewNodes = nc;
  for(i=0; i<nc; i++)
    mapMsg->newNodes[i] = computeMap->newNode(i);

  CProxy_WorkDistrib(thisgroup).recvComputeMapChanges(mapMsg);

/*
    // store the latest compute map
  SimParameters *simParams = Node::Object()->simParameters;
  if (simParams->storeComputeMap) {
    computeMap->saveComputeMap(simParams->computeMapFilename);
    CkPrintf("ComputeMap has been stored in %s.\n", simParams->computeMapFilename);
  }
*/
}

void WorkDistrib::recvComputeMapChanges(ComputeMapChangeMsg *msg) {
  
  if ( ! CkMyRank() ) {
    ComputeMap *computeMap = ComputeMap::Object();
    int i;
    for(i=0; i<computeMap->numComputes(); i++)
      computeMap->setNewNode(i,msg->newNodes[i]);
  }

  delete msg;

  CProxy_WorkDistrib workProxy(thisgroup);
  workProxy[0].doneSaveComputeMap();
}

void WorkDistrib::doneSaveComputeMap() {
  if (!--saveComputeMapCount) { 
    CkSendMsgBranch(saveComputeMapReturnEP, CkAllocMsg(0,0,0), 0, saveComputeMapReturnChareID);
  }
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////^M
//// Osman Sarood^M
//// Parallel Input Change^M
void WorkDistrib::setMapsArrived(bool s) {mapsArrived=s;}


#ifdef MEM_OPT_VERSION
void WorkDistrib::fillOnePatchCreationParallelIO(int patchId, FullAtomList *onePatchAtoms, Vector *velocities){
    

    CProxy_Node nd(CpvAccess(BOCclass_group).node);
    Node *node = nd.ckLocalBranch();
    PDB *pdb = node->pdb;
    SimParameters *params = node->simParameters;
    Molecule *molecule = node->molecule;

    const Lattice lattice = params->lattice;

    CProxy_PatchMgr pm(CpvAccess(BOCclass_group).patchMgr);
    PatchMgr *patchMgr = pm.ckLocalBranch();
    
    PatchMap *patchMap = PatchMap::Object();

    vector<int> *eachPatchAtomsList = Node::Object()->ioMgr->patchAtomList;
    vector<int> *thisPatchAtomsList = &eachPatchAtomsList[patchId];

    vector<int> *eachPatchAtomsList1 = Node::Object()->ioMgr->patchAtomIdMappingList;
    vector<int> *thisPatchAtomsList1 = &eachPatchAtomsList1[patchId];
    for(int i=0; i<thisPatchAtomsList->size(); i++){
        int aid = thisPatchAtomsList->at(i);
        int aidIdx=thisPatchAtomsList1->at(i);

        FullAtom a;
        a.id = aid;
                a.aidIdx=aidIdx;
                
        Position pos,vel;
        Node::Object()->ioMgr->get_position_for_patch_creation(&pos, aid,aidIdx);
                if ( Node::Object()->simParameters->initialTemp < 0.0 ) {
                        Node::Object()->ioMgr->get_position_for_patch_creationVel(&vel, aid,aidIdx);
                        a.velocity = vel;
                }
                else
                        a.velocity = velocities[aidIdx];        

        a.position = pos;
                

        
                

        #ifdef MEM_OPT_VERSION
/*
 *         a.sigId = molecule->getAtomSigId(aid);
 *                 a.exclId = molecule->getAtomExclSigId(aid);
 *                         a.vdwType = molecule->atomvdwtype(aid);
 *                         */
	FullAtom *at=Node::Object()->ioMgr->finalAtomList.begin();
/*
 *         a.sigId = Node::Object()->ioMgr->eachAtomSigPerProc[aidIdx];
 *                 a.exclId = Node::Object()->ioMgr->eachAtomExclSigPerProc[aidIdx];
 *                         a.vdwType = Node::Object()->ioMgr->vdwTypePerProc[aidIdx];
 *                         */
        a.sigId = at[aidIdx].sigId;
        a.exclId = at[aidIdx].exclId;
        a.vdwType = at[aidIdx].vdwType;

        #endif
        onePatchAtoms->add(a);

    }


    ScaledPosition center(0.5*(patchMap->min_a(patchId)+patchMap->max_a(patchId)),
                          0.5*(patchMap->min_b(patchId)+patchMap->max_b(patchId)),
                          0.5*(patchMap->min_c(patchId)+patchMap->max_c(patchId)));

    int n = onePatchAtoms->size();

    FullAtom *a = onePatchAtoms->begin();
    int j;

    Bool alchFepOn = params->alchFepOn;
    Bool alchThermIntOn = params->alchThermIntOn;


    Bool lesOn = params->lesOn;

    Bool pairInteractionOn = params->pairInteractionOn;

    Bool pressureProfileTypes = (params->pressureProfileAtomTypes > 1);
    Transform mother_transform;
    for(j=0; j < n; j++)
    {
      int aid = a[j].id;
                int aidIdx=a[j].aidIdx;
      if (params->splitPatch == SPLIT_PATCH_HYDROGEN) {
                if ( Node::Object()->ioMgr->is_hydrogenGroupParentPar(aidIdx) ) {
                        a[j].hydrogenGroupSize = (Node::Object()->ioMgr->is_hydrogenGroupParentPar(aidIdx)) ? Node::Object()->ioMgr->get_groupSizePar(aidIdx) : 0;
			a[j].migrationGroupSize= (Node::Object()->ioMgr->getIsMP(aidIdx)) ? Node::Object()->ioMgr->get_MigrationGpSizePar(aidIdx) : 0;
        } else {
          a[j].hydrogenGroupSize = 0;
	  a[j].migrationGroupSize=0;
        }
      } else {
        a[j].hydrogenGroupSize = 1;
      }


      a[j].nonbondedGroupSize = 0;  // must be set based on coordinates

      a[j].atomFixed = molecule->is_atom_fixed(aid) ? 1 : 0;
      a[j].fixedPosition = a[j].position;

      if ( a[j].hydrogenGroupSize ) {

        a[j].position = lattice.nearest(
                a[j].position, center, &(a[j].transform));

        mother_transform = a[j].transform;
      } else {
        a[j].position = lattice.apply_transform(a[j].position,mother_transform);
        a[j].transform = mother_transform;
      }
          Real *massPool=molecule->getAtomMassPool();
          Real *chargePool=molecule->getAtomChargePool();

      FullAtom *at=Node::Object()->ioMgr->finalAtomList.begin();
      a[j].mass = massPool[(Index)at[aidIdx].mass];

      a[j].charge = chargePool[(Index)at[aidIdx].charge];



      if ( alchFepOn || alchThermIntOn || lesOn || pairInteractionOn || pressureProfileTypes) {
        a[j].partition = molecule->get_fep_type(aid);
      }
      else {
        a[j].partition = 0;
      }

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
    if(params->fixedAtomsOn){
        int fixedCnt=0;
        for(j=0; j<n; j++)
            fixedCnt += (a[j].atomFixed ? 1:0);
        patchMgr->setHomePatchFixedAtomNum(patchId, fixedCnt);
    }

    if ( params->outputPatchDetails ) {
      int numAtomsInPatch = n;
      int numFixedAtomsInPatch = 0;
      int numAtomsInFixedGroupsInPatch = 0;
      for(j=0; j < n; j++) {
        numFixedAtomsInPatch += ( a[j].atomFixed ? 1 : 0 );
        numAtomsInFixedGroupsInPatch += ( a[j].groupFixed ? 1 : 0 );
      }
      iout << "PATCH_DETAILS:"
           << " patch " << patchId
           << " atoms " << numAtomsInPatch
           << " fixed_atoms " << numFixedAtomsInPatch
           << " fixed_groups " << numAtomsInFixedGroupsInPatch
           << "\n" << endi;
    }
}


void WorkDistrib::random_velocities_parallel(BigReal Temp,Molecule *structure,
                                    Vector *v,int newTotalAtoms)
{
  int i, j;             //  Loop counter
  BigReal kbT;          //  Boltzman constant * Temp
  BigReal randnum;      //  Random number from -6.0 to 6.0
  BigReal kbToverM;     //  sqrt(Kb*Temp/Mass)
  SimParameters *simParams = Node::Object()->simParameters;
  Bool lesOn = simParams->lesOn;
  Random vel_random(simParams->randomSeed);
  PDB *pdb=Node::Object()->pdb;
  int lesReduceTemp = lesOn && simParams->lesReduceTemp;
  BigReal tempFactor = lesReduceTemp ? 1.0 / simParams->lesFactor : 1.0;

  kbT = Temp*BOLTZMAN;
  int count=0;
  for(int i=0;i<newTotalAtoms;i++)
  {
    FullAtom *at=Node::Object()->ioMgr->finalAtomList.begin();
    int aid=at[i].id;

    Real *massPool=Node::Object()->molecule->getAtomMassPool();
    Real atomMs=massPool[(Index)at[i].mass];

    if (atomMs <= 0.) {
      kbToverM = 0.;
    } else {
      kbToverM = sqrt(kbT *
        ( lesOn && structure->get_fep_type(aid) ? tempFactor : 1.0 ) /
                          atomMs );
    }
    for (randnum=0.0, j=0; j<12; j++)
    {
      randnum += vel_random.uniform();
    }

    randnum -= 6.0;

    v[count].x = randnum*kbToverM;

    for (randnum=0.0, j=0; j<12; j++)
    {
      randnum += vel_random.uniform();
    }

    randnum -= 6.0;

    v[count].y = randnum*kbToverM;

    for (randnum=0.0, j=0; j<12; j++)
    {
      randnum += vel_random.uniform();
    }

    randnum -= 6.0;
    
    v[count].z = randnum*kbToverM;
    count++;
  }


}
#endif
/////////////////////////////////////////////////////////////////////////////////////////////////////////////

//only called on node 0
int *WorkDistrib::caclNumAtomsInEachPatch(){
    StringList *current;
    int i;
    CProxy_Node nd(CkpvAccess(BOCclass_group).node);
    Node *node = nd.ckLocalBranch();
    PatchMap *patchMap = PatchMap::Object();
    CProxy_PatchMgr pm(CkpvAccess(BOCclass_group).patchMgr);
    PatchMgr *patchMgr = pm.ckLocalBranch();
    SimParameters *params = node->simParameters;
    Molecule *molecule = node->molecule;
    PDB *pdb = node->pdb;

    int numPatches = patchMap->numPatches();
    int numAtoms = pdb->num_atoms();    

    vector<int> *eachPatchAtomList = patchMap->getTmpPatchAtomsList();

    int *patchAtomCnt = new int[numPatches];    
    memset(patchAtomCnt, 0, sizeof(int)*numPatches);

    const Lattice lattice = params->lattice;

    Position eachAtomPos;
      {
      // split atoms into patched based on helix-group and position
      int aid, pid=0;
      for(i=0; i < numAtoms; i++)
        {        
        // Assign atoms to patches without splitting hydrogen groups.
        // We know that the hydrogenGroup array is sorted with group parents
        // listed first.  Thus, only change the pid if an atom is a group parent.
        aid = molecule->hydrogenGroup[i].atomID;
        pdb->get_position_for_atom(&eachAtomPos, aid);
        if (molecule->hydrogenGroup[i].isMP)            
            pid = patchMap->assignToPatch(eachAtomPos,lattice);
        // else: don't change pid        
        patchAtomCnt[pid]++;
        eachPatchAtomList[pid].push_back(aid);
        }
      }

    return patchAtomCnt;    
}

//----------------------------------------------------------------------
// This should only be called on node 0.
//----------------------------------------------------------------------
FullAtomList *WorkDistrib::createAtomLists(void)
{
  int i;
  StringList *current;	//  Pointer used to retrieve configuration items
  CProxy_Node nd(CkpvAccess(BOCclass_group).node);
  Node *node = nd.ckLocalBranch();
  PatchMap *patchMap = PatchMap::Object();
  CProxy_PatchMgr pm(CkpvAccess(BOCclass_group).patchMgr);
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

  //  If COMMotion == no, remove center of mass motion
  if (!(params->comMove)) {
    remove_com_motion(velocities, molecule, numAtoms);
  }

  FullAtomList *atoms = new FullAtomList[numPatches];

  const Lattice lattice = params->lattice;

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
      FullAtom a;
      a.id = aid;
      a.position = positions[aid];
      a.velocity = velocities[aid];
      #ifdef MEM_OPT_VERSION
      a.sigId = molecule->getAtomSigId(aid);
      a.exclId = molecule->getAtomExclSigId(aid);
      #endif
      a.vdwType = molecule->atomvdwtype(aid);
      HydrogenGroupID &h = molecule->hydrogenGroup[i];
      a.hydrogenGroupSize = h.isGP ? h.atomsInGroup : 0;
      a.migrationGroupSize = h.isMP ? h.atomsInMigrationGroup : 0;
      if (h.isMP) {
	pid = patchMap->assignToPatch(positions[aid],lattice);
      } // else: don't change pid
      atoms[pid].add(a);
      }
    }

  delete [] positions;
  delete [] velocities;

  for(i=0; i < numPatches; i++)
  {
    ScaledPosition center(0.5*(patchMap->min_a(i)+patchMap->max_a(i)),
			  0.5*(patchMap->min_b(i)+patchMap->max_b(i)),
			  0.5*(patchMap->min_c(i)+patchMap->max_c(i)));

    int n = atoms[i].size();
    FullAtom *a = atoms[i].begin();
    int j;
//Modifications for alchemical fep
    Bool alchFepOn = params->alchFepOn;
    Bool alchThermIntOn = params->alchThermIntOn;
//fepe
    Bool lesOn = params->lesOn;
  
    Bool pairInteractionOn = params->pairInteractionOn;

    Bool pressureProfileTypes = (params->pressureProfileAtomTypes > 1);

    Transform mother_transform;
    for(j=0; j < n; j++)
    {
      int aid = a[j].id;

      a[j].nonbondedGroupSize = 0;  // must be set based on coordinates

      a[j].atomFixed = molecule->is_atom_fixed(aid) ? 1 : 0;
      a[j].fixedPosition = a[j].position;

      if ( a[j].hydrogenGroupSize ) {
        a[j].position = lattice.nearest(
		a[j].position, center, &(a[j].transform));
        mother_transform = a[j].transform;
      } else {
        a[j].position = lattice.apply_transform(a[j].position,mother_transform);
        a[j].transform = mother_transform;
      }

      a[j].mass = molecule->atommass(aid);
      a[j].charge = molecule->atomcharge(aid);

//Modifications for alchemical fep
      if ( alchFepOn || alchThermIntOn || lesOn || pairInteractionOn || pressureProfileTypes) {
        a[j].partition = molecule->get_fep_type(aid);
      } 
      else {
        a[j].partition = 0;
      }
//fepe

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

    if ( params->outputPatchDetails ) {
      int patchId = i;
      int numAtomsInPatch = n;
      int numFixedAtomsInPatch = 0;
      int numAtomsInFixedGroupsInPatch = 0;
      for(j=0; j < n; j++) {
        numFixedAtomsInPatch += ( a[j].atomFixed ? 1 : 0 );
        numAtomsInFixedGroupsInPatch += ( a[j].groupFixed ? 1 : 0 );
      }
      iout << "PATCH_DETAILS:"
           << " patch " << patchId
           << " atoms " << numAtomsInPatch
           << " fixed_atoms " << numFixedAtomsInPatch
           << " fixed_groups " << numAtomsInFixedGroupsInPatch
           << "\n" << endi;
    }
  }

  return atoms;

}

//This should only be called on node 0
//This differs from creatHomePatches in: create home patches
//without populating them with actual atoms' data. The only field
//set is the number of atoms each patch contains.
void WorkDistrib::preCreateHomePatches(){

  PatchMap *patchMap = PatchMap::Object();
  CProxy_PatchMgr pm(CkpvAccess(BOCclass_group).patchMgr);
  PatchMgr *patchMgr = pm.ckLocalBranch();

  int numPatches = patchMap->numPatches();

  patchMap->initTmpPatchAtomsList();

  int *patchAtomCnt = caclNumAtomsInEachPatch();

  int maxAtoms = -1;
  int maxPatch = -1;
  for(int i=0; i < numPatches; i++) {
    int numAtoms = patchAtomCnt[i];
    if ( numAtoms > maxAtoms ) { maxAtoms = numAtoms; maxPatch = i; }
  }
  iout << iINFO << "LARGEST PATCH (" << maxPatch <<
	") HAS " << maxAtoms << " ATOMS\n" << endi;

  for(int i=0; i < numPatches; i++)
  {
    if ( ! ( i % 100 ) )
    {
      DebugM(3,"Pre-created " << i << " patches so far.\n");
    }

    patchMgr->preCreateHomePatch(i,patchAtomCnt[i]);
  }

  delete [] patchAtomCnt;
}


//----------------------------------------------------------------------
// This should only be called on node 0.
//----------------------------------------------------------------------
void WorkDistrib::createHomePatches(void)
{
  int i;
  PatchMap *patchMap = PatchMap::Object();
  CProxy_PatchMgr pm(CkpvAccess(BOCclass_group).patchMgr);
  PatchMgr *patchMgr = pm.ckLocalBranch();

  int numPatches = patchMap->numPatches();

  FullAtomList *atoms = createAtomLists();
    
#ifdef MEM_OPT_VERSION
/*  CProxy_Node nd(CkpvAccess(BOCclass_group).node);
  Node *node = nd.ckLocalBranch();
  node->molecule->delEachAtomSigs();
  node->molecule->delMassChargeSpace();
*/
#endif

  int maxAtoms = -1;
  int maxPatch = -1;
  for(i=0; i < numPatches; i++) {
    int numAtoms = atoms[i].size();
    if ( numAtoms > maxAtoms ) { maxAtoms = numAtoms; maxPatch = i; }
  }
  iout << iINFO << "LARGEST PATCH (" << maxPatch <<
	") HAS " << maxAtoms << " ATOMS\n" << endi;

  for(i=0; i < numPatches; i++)
  {
    if ( ! ( i % 100 ) )
    {
      DebugM(3,"Created " << i << " patches so far.\n");
    }

    patchMgr->createHomePatch(i,atoms[i]);
  }

  delete [] atoms;
}

void WorkDistrib::fillOnePatchAtoms(int patchId, FullAtomList *onePatchAtoms, Vector *velocities){
    // ref BOC
    CProxy_Node nd(CkpvAccess(BOCclass_group).node);
    Node *node = nd.ckLocalBranch();
    PDB *pdb = node->pdb;
    SimParameters *params = node->simParameters;
    Molecule *molecule = node->molecule;

    const Lattice lattice = params->lattice;

    CProxy_PatchMgr pm(CkpvAccess(BOCclass_group).patchMgr);
    PatchMgr *patchMgr = pm.ckLocalBranch();
    // ref singleton
    PatchMap *patchMap = PatchMap::Object();

    vector<int> *eachPatchAtomsList = patchMap->getTmpPatchAtomsList();
    vector<int> *thisPatchAtomsList = &eachPatchAtomsList[patchId];

    for(int i=0; i<thisPatchAtomsList->size(); i++){
        int aid = thisPatchAtomsList->at(i);
        FullAtom a;
        a.id = aid;
        Position pos;
        pdb->get_position_for_atom(&pos, aid);
        a.position = pos;
        a.velocity = velocities[aid];
        #ifdef MEM_OPT_VERSION
        a.sigId = molecule->getAtomSigId(aid);
        a.exclId = molecule->getAtomExclSigId(aid);
        #endif
        a.vdwType = molecule->atomvdwtype(aid);
        HydrogenGroupID &h =
          molecule->hydrogenGroup[molecule->atoms[aid].hydrogenList];
        a.hydrogenGroupSize = h.isGP ? h.atomsInGroup : 0;
        a.migrationGroupSize = h.isMP ? h.atomsInMigrationGroup : 0;
        onePatchAtoms->add(a);
    }

    ScaledPosition center(0.5*(patchMap->min_a(patchId)+patchMap->max_a(patchId)),
                          0.5*(patchMap->min_b(patchId)+patchMap->max_b(patchId)),
                          0.5*(patchMap->min_c(patchId)+patchMap->max_c(patchId)));
    int n = onePatchAtoms->size();
    FullAtom *a = onePatchAtoms->begin();
    int j;
//Modifications for alchemical fep
    Bool alchFepOn = params->alchFepOn;
    Bool alchThermIntOn = params->alchThermIntOn;
//fepe
    Bool lesOn = params->lesOn;
  
    Bool pairInteractionOn = params->pairInteractionOn;

    Bool pressureProfileTypes = (params->pressureProfileAtomTypes > 1);

    Transform mother_transform;
    for(j=0; j < n; j++)
    {
      int aid = a[j].id;

      a[j].nonbondedGroupSize = 0;  // must be set based on coordinates

      a[j].atomFixed = molecule->is_atom_fixed(aid) ? 1 : 0;
      a[j].fixedPosition = a[j].position;

      if ( a[j].hydrogenGroupSize ) {
        a[j].position = lattice.nearest(
		a[j].position, center, &(a[j].transform));
        mother_transform = a[j].transform;
      } else {
        a[j].position = lattice.apply_transform(a[j].position,mother_transform);
        a[j].transform = mother_transform;
      }

      a[j].mass = molecule->atommass(aid);
      a[j].charge = molecule->atomcharge(aid);

//Modifications for alchemical fep
      if ( alchFepOn || alchThermIntOn || lesOn || pairInteractionOn || pressureProfileTypes) {
        a[j].partition = molecule->get_fep_type(aid);
      } 
      else {
        a[j].partition = 0;
      }
//fepe
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

    if(params->fixedAtomsOn){
	int fixedCnt=0;
	for(j=0; j<n; j++)
	    fixedCnt += (a[j].atomFixed ? 1:0);
	patchMgr->setHomePatchFixedAtomNum(patchId, fixedCnt);
    }    

    if ( params->outputPatchDetails ) {    
      int numAtomsInPatch = n;
      int numFixedAtomsInPatch = 0;
      int numAtomsInFixedGroupsInPatch = 0;
      for(j=0; j < n; j++) {
        numFixedAtomsInPatch += ( a[j].atomFixed ? 1 : 0 );
        numAtomsInFixedGroupsInPatch += ( a[j].groupFixed ? 1 : 0 );
      }
      iout << "PATCH_DETAILS:"
           << " patch " << patchId
           << " atoms " << numAtomsInPatch
           << " fixed_atoms " << numFixedAtomsInPatch
           << " fixed_groups " << numAtomsInFixedGroupsInPatch
           << "\n" << endi;
    }
}

//should be called only on node 0
void WorkDistrib::initAndSendHomePatch(){
  StringList *current;
  // ref BOC
  CProxy_Node nd(CkpvAccess(BOCclass_group).node);
  Node *node = nd.ckLocalBranch();
  Molecule *molecule = node->molecule;
  PDB *pdb = node->pdb;
  SimParameters *params = node->simParameters;

  CProxy_PatchMgr pm(CkpvAccess(BOCclass_group).patchMgr);
  PatchMgr *patchMgr = pm.ckLocalBranch();
  // ref singleton
  PatchMap *patchMap = PatchMap::Object();

  int numAtoms = pdb->num_atoms();

  //1. create atoms' velocities
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

  //  If COMMotion == no, remove center of mass motion
  if (!(params->comMove)) {
    remove_com_motion(velocities, molecule, numAtoms);
  }


  for(int i=0; i<patchMap->numPatches(); i++){
      //2. fill each patch with actual atom data
      //the work flow should looks like the createAtomList
      FullAtomList *onePatchAtoms = new FullAtomList;
      fillOnePatchAtoms(i, onePatchAtoms, velocities);
      
      patchMgr->fillHomePatchAtomList(i, onePatchAtoms);

      delete onePatchAtoms;

      //distribute this home patch
      if(patchMap->node(i) != node->myid()){
          //need to move to other nodes
          patchMgr->sendOneHomePatch(i, patchMap->node(i));
      }
      
  }  

  delete [] velocities;
  patchMap->delTmpPatchAtomsList();
}

void WorkDistrib::distributeHomePatches() {
  // ref BOC
  CProxy_Node nd(CkpvAccess(BOCclass_group).node);
  Node *node = nd.ckLocalBranch();
  CProxy_PatchMgr pm(CkpvAccess(BOCclass_group).patchMgr);
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

void WorkDistrib::reinitAtoms() {

  PatchMap *patchMap = PatchMap::Object();
  CProxy_PatchMgr pm(CkpvAccess(BOCclass_group).patchMgr);
  PatchMgr *patchMgr = pm.ckLocalBranch();

  int numPatches = patchMap->numPatches();

  FullAtomList *atoms = createAtomLists();

  for(int i=0; i < numPatches; i++) {
    patchMgr->sendAtoms(i,atoms[i]);
  }

  delete [] atoms;

}


//----------------------------------------------------------------------

class MapDistribMsg: public CMessage_MapDistribMsg {
  public:
    char *patchMapData;
    char *computeMapData;

//  VARSIZE_DECL(MapDistribMsg);
};

/*
VARSIZE_MSG(MapDistribMsg,
  VARSIZE_ARRAY(patchMapData);
  VARSIZE_ARRAY(computeMapData);
)
*/

void WorkDistrib::sendMaps(void)
{
  if ( CkNumPes() == 1 ) {
    mapsArrived = true;
    return;
  }

  //Automatically enable spanning tree
  CProxy_Node nd(CkpvAccess(BOCclass_group).node);
  Node *node = nd.ckLocalBranch();
  SimParameters *params = node->simParameters;
  if(PatchMap::Object()->numPatches() <= CkNumPes()/4 && params->isSendSpanningTreeUnset())
    ProxyMgr::Object()->setSendSpanning();

  int sizes[2];
  sizes[0] = PatchMap::Object()->packSize();
  sizes[1] = ComputeMap::Object()->packSize();

  MapDistribMsg *mapMsg = new (sizes[0], sizes[1], 0) MapDistribMsg;

  PatchMap::Object()->pack(mapMsg->patchMapData);
  ComputeMap::Object()->pack(mapMsg->computeMapData);

  CProxy_WorkDistrib workProxy(thisgroup);
  workProxy[0].saveMaps(mapMsg);
}

// saveMaps() is called when the map message is received
void WorkDistrib::saveMaps(MapDistribMsg *msg)
{
  // Use a resend to forward messages before processing.  Otherwise the
  // map distribution is slow on many CPUs.  We need to use a tree
  // rather than a broadcast because some implementations of broadcast
  // generate a copy of the message on the sender for each recipient.
  // This is because MPI doesn't allow re-use of an outstanding buffer.

  if ( mapsArrived && CkMyPe() ) {
    PatchMap::Object()->unpack(msg->patchMapData);
    if ( ! CkMyRank() ) ComputeMap::Object()->unpack(msg->computeMapData);

    //Automatically enable spanning tree
    CProxy_Node nd(CkpvAccess(BOCclass_group).node);
    Node *node = nd.ckLocalBranch();
    SimParameters *params = node->simParameters;
    if(PatchMap::Object()->numPatches() <= CkNumPes()/4 && params->isSendSpanningTreeUnset())
      ProxyMgr::Object()->setSendSpanning();
  }
  if ( mapsArrived ) {
    delete msg;
    return;
  }

  mapsArrived = true;

  int pids[3];
  int basePe = 2 * CkMyPe() + 1;
  int npid = 0;
  if ( (basePe+npid) < CkNumPes() ) { pids[npid] = basePe + npid; ++npid; }
  if ( (basePe+npid) < CkNumPes() ) { pids[npid] = basePe + npid; ++npid; }
  pids[npid] = CkMyPe(); ++npid;  // always send the message to ourselves
  CProxy_WorkDistrib(thisgroup).saveMaps(msg,npid,pids);
}


//----------------------------------------------------------------------
void WorkDistrib::patchMapInit(void)
{
  PatchMap *patchMap = PatchMap::Object();
  CProxy_Node nd(CkpvAccess(BOCclass_group).node);
  Node *node = nd.ckLocalBranch();
  SimParameters *params = node->simParameters;
  Lattice lattice = params->lattice;

  BigReal patchSize = params->patchDimension;

  ScaledPosition xmin, xmax;

  double maxNumPatches = 1000000;  // need to adjust fractional values
  if ( params->minAtomsPerPatch > 0 )
#ifndef MEM_OPT_VERSION
    maxNumPatches = node->pdb->num_atoms() / params->minAtomsPerPatch;
#else
    maxNumPatches = Node::Object()->ioMgr->totalAtomsSys / params->minAtomsPerPatch;
#endif

  DebugM(3,"Mapping patches\n");
  if ( lattice.a_p() && lattice.b_p() && lattice.c_p() ) {
    xmin = 0.;  xmax = 0.;
  } else {
  // Need to use full box for FMA to match NAMD 1.X results.
  if ( params->FMAOn ) {
    node->pdb->find_extremes(&(xmin.x),&(xmax.x),lattice.a_r());
    node->pdb->find_extremes(&(xmin.y),&(xmax.y),lattice.b_r());
    node->pdb->find_extremes(&(xmin.z),&(xmax.z),lattice.c_r());
  // Otherwise, this allows a small number of stray atoms.
  } else {
    node->pdb->find_extremes(&(xmin.x),&(xmax.x),lattice.a_r(),0.9);
    node->pdb->find_extremes(&(xmin.y),&(xmax.y),lattice.b_r(),0.9);
    node->pdb->find_extremes(&(xmin.z),&(xmax.z),lattice.c_r(),0.9);
  }
  }

  BigReal origin_shift;
  origin_shift = lattice.a_r() * lattice.origin();
  xmin.x -= origin_shift;
  xmax.x -= origin_shift;
  origin_shift = lattice.b_r() * lattice.origin();
  xmin.y -= origin_shift;
  xmax.y -= origin_shift;
  origin_shift = lattice.c_r() * lattice.origin();
  xmin.z -= origin_shift;
  xmax.z -= origin_shift;

  // SimParameters default is -1 for unset
  int twoAwayX = params->twoAwayX;
  int twoAwayY = params->twoAwayY;
  int twoAwayZ = params->twoAwayZ;

#ifdef NAMD_CUDA
  // for CUDA be sure there are more patches than pes

  int numPatches = patchMap->sizeGrid(
	xmin,xmax,lattice,patchSize,maxNumPatches,
	twoAwayX>0 ? 2 : 1, twoAwayY>0 ? 2 : 1, twoAwayZ>0 ? 2 : 1);
  if ( numPatches < CkNumPes() && twoAwayX < 0 ) {
    twoAwayX = 1;
    numPatches = patchMap->sizeGrid(
	xmin,xmax,lattice,patchSize,maxNumPatches,
	twoAwayX>0 ? 2 : 1, twoAwayY>0 ? 2 : 1, twoAwayZ>0 ? 2 : 1);
  }
  if ( numPatches < CkNumPes() && twoAwayY < 0 ) {
    twoAwayY = 1;
    numPatches = patchMap->sizeGrid(
	xmin,xmax,lattice,patchSize,maxNumPatches,
	twoAwayX>0 ? 2 : 1, twoAwayY>0 ? 2 : 1, twoAwayZ>0 ? 2 : 1);
  }
  if ( numPatches < CkNumPes() && twoAwayZ < 0 ) {
    twoAwayZ = 1;
    numPatches = patchMap->sizeGrid(
	xmin,xmax,lattice,patchSize,maxNumPatches,
	twoAwayX>0 ? 2 : 1, twoAwayY>0 ? 2 : 1, twoAwayZ>0 ? 2 : 1);
  }
  if ( numPatches < CkNumPes() ) {
    NAMD_die("CUDA-enabled NAMD requires more patches than processes.");
  }

  patchMap->makePatches(xmin,xmax,lattice,patchSize,maxNumPatches,
	twoAwayX>0 ? 2 : 1, twoAwayY>0 ? 2 : 1, twoAwayZ>0 ? 2 : 1);

#else

  int numPatches = patchMap->sizeGrid(
	xmin,xmax,lattice,patchSize,maxNumPatches,
	twoAwayX ? 2 : 1, twoAwayY ? 2 : 1, twoAwayZ ? 2 : 1);
  if ( numPatches > (CkNumPes() - 1) && twoAwayZ < 0 ) {
    twoAwayZ = 0;
    numPatches = patchMap->sizeGrid(
	xmin,xmax,lattice,patchSize,maxNumPatches,
	twoAwayX ? 2 : 1, twoAwayY ? 2 : 1, twoAwayZ ? 2 : 1);
  }
  if ( numPatches > (CkNumPes() - 1) && twoAwayY < 0 ) {
    twoAwayY = 0;
    numPatches = patchMap->sizeGrid(
	xmin,xmax,lattice,patchSize,maxNumPatches,
	twoAwayX ? 2 : 1, twoAwayY ? 2 : 1, twoAwayZ ? 2 : 1);
  }
  if ( numPatches > (CkNumPes() - 1) && twoAwayX < 0 ) {
    twoAwayX = 0;
    numPatches = patchMap->sizeGrid(
	xmin,xmax,lattice,patchSize,maxNumPatches,
	twoAwayX ? 2 : 1, twoAwayY ? 2 : 1, twoAwayZ ? 2 : 1);
  }

  patchMap->makePatches(xmin,xmax,lattice,patchSize,maxNumPatches,
	twoAwayX ? 2 : 1, twoAwayY ? 2 : 1, twoAwayZ ? 2 : 1);

#endif

}


//----------------------------------------------------------------------
void WorkDistrib::assignNodeToPatch()
{
  PatchMap *patchMap = PatchMap::Object();
  int nNodes = Node::Object()->numNodes();

#if USE_TOPOMAP 
  TopoManager tmgr;
  int numPes = tmgr.getDimNX() * tmgr.getDimNY() * tmgr.getDimNZ();
  if (numPes > patchMap->numPatches() && (assignPatchesTopoGridRecBisection() > 0)) {
    CkPrintf ("Blue Gene/L topology partitioner finished successfully \n");
  }
  else
#endif
    if (nNodes > patchMap->numPatches())
      assignPatchesBitReversal();
    else
#if (CMK_BLUEGENEP | CMK_BLUEGENEL)
      assignPatchesRecursiveBisection();
#else
      assignPatchesSpaceFillingCurve();
#endif
      // assignPatchesRoundRobin();
      // assignPatchesToLowestLoadNode();
  
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
printf("numAtoms=%d Node::Object()->molecule->numAtoms=%d\n",numAtoms,Node::Object()->molecule->numAtoms);
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
//   Node *node = CLocalBranch(Node, CkpvAccess(BOCclass_group).node);

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
  CProxy_Node nd(CkpvAccess(BOCclass_group).node);
  Node *node = nd.ckLocalBranch();
  SimParameters *simParams = node->simParameters;

  int *load = new int[node->numNodes()];
  int *assignedNodes = new int[patchMap->numPatches()];
  for (int i=0; i<node->numNodes(); i++) {
    load[i] = 0;
  }

  int defaultNode = 0;
  if ( simParams->noPatchesOnZero && node->numNodes() > 1 ){
    defaultNode = 1;
    if( simParams->noPatchesOnOne && node->numNodes() > 2)
      defaultNode = 2;
  }
  // Assign patch to node with least atoms assigned.
  for(pid=0; pid < patchMap->numPatches(); pid++) {
    assignedNode = defaultNode;
    for (int i=assignedNode + 1; i < node->numNodes(); i++) {
      if (load[i] < load[assignedNode]) assignedNode = i;
    }
    assignedNodes[pid] = assignedNode;
    load[assignedNode] += patchMap->patch(pid)->getNumAtoms() + 1;
  }

  delete[] load;
  sortNodesAndAssign(assignedNodes);
  delete[] assignedNodes;
}

//----------------------------------------------------------------------
void WorkDistrib::assignPatchesBitReversal() 
{
  int pid; 
  PatchMap *patchMap = PatchMap::Object();
  CProxy_Node nd(CkpvAccess(BOCclass_group).node);
  Node *node = nd.ckLocalBranch();

  int ncpus = node->numNodes();
  int npatches = patchMap->numPatches();
  if ( ncpus <= npatches )
    NAMD_bug("WorkDistrib::assignPatchesBitReversal called improperly");

  // find next highest power of two
  int npow2 = 1;  int nbits = 0;
  while ( npow2 < ncpus ) { npow2 *= 2; nbits += 1; }

  // build bit reversal sequence
  SortableResizeArray<int> seq(ncpus);
  // avoid using node 0 (reverse of 0 is 0 so start at 1)
  int i = 1;
  for ( int icpu=0; icpu<(ncpus-1); ++icpu ) {
    int ri;
    for ( ri = ncpus; ri >= ncpus; ++i ) {
      ri = 0;
      int pow2 = 1;
      int rpow2 = npow2 / 2;
      for ( int j=0; j<nbits; ++j ) {
        ri += rpow2 * ( ( i / pow2 ) % 2 );
        pow2 *= 2;  rpow2 /= 2;
      }
    }
    seq[icpu] = ri;
  }

  // extract and sort patch locations
  sortNodesAndAssign(seq.begin());
  if ( ncpus > 2*npatches ) sortNodesAndAssign(seq.begin()+npatches, 1);
}

//----------------------------------------------------------------------
struct nodesort {
  int node;
  int a_total;
  int b_total;
  int c_total;
  int npatches;
  nodesort() : node(-1),a_total(0),b_total(0),c_total(0),npatches(0) { ; }
  int operator==(const nodesort &o) const {
    float a1 = ((float)a_total)/((float)npatches);
    float a2 = ((float)o.a_total)/((float)o.npatches);
    float b1 = ((float)b_total)/((float)npatches);
    float b2 = ((float)o.b_total)/((float)o.npatches);
    float c1 = ((float)c_total)/((float)npatches);
    float c2 = ((float)o.c_total)/((float)o.npatches);
    return ((a1 == a2) && (b1 == b2) && (c1 == c2));
  }
  int operator<(const nodesort &o) const {
    float a1 = ((float)a_total)/((float)npatches);
    float a2 = ((float)o.a_total)/((float)o.npatches);
    float b1 = ((float)b_total)/((float)npatches);
    float b2 = ((float)o.b_total)/((float)o.npatches);
    float c1 = ((float)c_total)/((float)npatches);
    float c2 = ((float)o.c_total)/((float)o.npatches);
    return ( (a1 < a2) || ((a1 == a2) && (b1 < b2)) ||
		((a1 == a2) && (b1 == b2) && (c1 < c2)) );
  }
};


void WorkDistrib::sortNodesAndAssign(int *assignedNode, int baseNodes) {
  // if baseNodes is zero (default) then set both nodes and basenodes
  // if baseNodes is nonzero then this is a second call to set basenodes only
  int i, pid; 
  PatchMap *patchMap = PatchMap::Object();
  CProxy_Node nd(CkpvAccess(BOCclass_group).node);
  Node *node = nd.ckLocalBranch();
  int nnodes = node->numNodes();
  int npatches = patchMap->numPatches();

  ResizeArray<nodesort> allnodes(nnodes);
  for ( i=0; i < nnodes; ++i ) {
    allnodes[i].node = i;
  }
  for ( pid=0; pid<npatches; ++pid ) {
    // iout << pid << " " << assignedNode[pid] << "\n" << endi;
    allnodes[assignedNode[pid]].npatches++;
    allnodes[assignedNode[pid]].a_total += patchMap->index_a(pid);
    allnodes[assignedNode[pid]].b_total += patchMap->index_b(pid);
    allnodes[assignedNode[pid]].c_total += patchMap->index_c(pid);
  }
  SortableResizeArray<nodesort> usednodes(nnodes);
  usednodes.resize(0);
  for ( i=0; i < nnodes; ++i ) {
    if ( allnodes[i].npatches ) usednodes.add(allnodes[i]);
  }
  usednodes.sort();
  int nused = usednodes.size();
  int i2 = nused/2;
  for ( i=0; i < nnodes; ++i ) {
    if ( allnodes[i].npatches ) allnodes[usednodes[i2++].node].node = i;
    if ( i2 == nused ) i2 = 0;
  }

  for ( pid=0; pid<npatches; ++pid ) {
    // iout << pid << " " <<  allnodes[assignedNode[pid]].node << "\n" << endi;
    if ( ! baseNodes ) {
      patchMap->assignNode(pid, allnodes[assignedNode[pid]].node);
    }
    patchMap->assignBaseNode(pid, allnodes[assignedNode[pid]].node);
  }
}

//----------------------------------------------------------------------
void WorkDistrib::assignPatchesRoundRobin() 
{
  int pid; 
  PatchMap *patchMap = PatchMap::Object();
  CProxy_Node nd(CkpvAccess(BOCclass_group).node);
  Node *node = nd.ckLocalBranch();
  int *assignedNode = new int[patchMap->numPatches()];

  for(pid=0; pid < patchMap->numPatches(); pid++) {
    assignedNode[pid] = pid % node->numNodes();
  }

  sortNodesAndAssign(assignedNode);
  delete [] assignedNode;
}

//----------------------------------------------------------------------
void WorkDistrib::assignPatchesRecursiveBisection() 
{
  PatchMap *patchMap = PatchMap::Object();
  int *assignedNode = new int[patchMap->numPatches()];
  int numNodes = Node::Object()->numNodes();
  SimParameters *simParams = Node::Object()->simParameters;
  int usedNodes = numNodes;
  int unusedNodes = 0;
  if ( simParams->noPatchesOnZero && numNodes > 1 ){
    usedNodes -= 1;
    if(simParams->noPatchesOnOne && numNodes > 2)
      usedNodes -= 1;
  }  
  unusedNodes = numNodes - usedNodes;
  RecBisection recBisec(usedNodes,PatchMap::Object());
  if ( recBisec.partition(assignedNode) ) {
    if ( unusedNodes !=0 ) {
      for ( int i=0; i<patchMap->numPatches(); ++i ) {
        assignedNode[i] += unusedNodes;
      }
    }
    sortNodesAndAssign(assignedNode);
    delete [] assignedNode;
  } else {
    //free the array here since a same array will be allocated
    //in assignPatchesToLowestLoadNode function, thus reducting
    //temporary memory usage
    delete [] assignedNode; 
    
    iout << iWARN 
	 << "WorkDistrib: Recursive bisection fails,"
	 << "invoking least-load algorithm\n";
    assignPatchesToLowestLoadNode();
  }
}

//----------------------------------------------------------------------
void WorkDistrib::assignPatchesSpaceFillingCurve() 
{
  PatchMap *patchMap = PatchMap::Object();
  int *assignedNode = new int[patchMap->numPatches()];
  int numNodes = Node::Object()->numNodes();
  SimParameters *simParams = Node::Object()->simParameters;
  int usedNodes = numNodes;
  int unusedNodes = 0;
  if ( simParams->noPatchesOnZero && numNodes > 1 ){
    usedNodes -= 1;
    if(simParams->noPatchesOnOne && numNodes > 2)
      usedNodes -= 1;
  }  
  unusedNodes = numNodes - usedNodes;

  int numPatches = patchMap->numPatches();
  if ( numPatches < usedNodes )
    NAMD_bug("WorkDistrib::assignPatchesSpaceFillingCurve() called with more nodes than patches");

  ResizeArray<double> patchLoads(numPatches);
  SortableResizeArray<double> sortedLoads(numPatches);
  for ( int i=0; i<numPatches; ++i ) {
    double load = patchMap->patch(i)->getNumAtoms() + 10;
    patchLoads[i] = load;
    sortedLoads[i] = load;
  }
  sortedLoads.sort();

  // limit maxPatchLoad to adjusted average load per node
  double sumLoad = 0;
  double maxPatchLoad = 1;
  for ( int i=0; i<numPatches; ++i ) {
    double load = sortedLoads[i];
    double total = sumLoad + (numPatches-i) * load;
    if ( usedNodes * load > total ) break;
    sumLoad += load;
    maxPatchLoad = load;
  }
  double totalLoad = 0;
  for ( int i=0; i<numPatches; ++i ) {
    if ( patchLoads[i] > maxPatchLoad ) patchLoads[i] = maxPatchLoad;
    totalLoad += patchLoads[i];
  }
  if ( usedNodes * maxPatchLoad > totalLoad )
    NAMD_bug("algorithm failure in WorkDistrib::assignPatchesSpaceFillingCurve()");

  // walk through patches in space-filling curve
  sumLoad = 0;
  int node = 0;
  int adim = patchMap->gridsize_a();
  int bdim = patchMap->gridsize_b();
  int cdim = patchMap->gridsize_c();
  int b = 0;
  int c = 0;
  int binc = 1;
  int cinc = 1;
  for ( int a = 0; a < adim; ++a ) {
    while ( b >= 0 && b < bdim ) {
      while ( c >= 0 && c < cdim ) {
        int pid = patchMap->pid(a,b,c);
        assignedNode[pid] = node;
        sumLoad += patchLoads[pid];
        double targetLoad = (double)(node+1) / (double)usedNodes;
        targetLoad *= totalLoad;
        if ( node+1 < usedNodes && sumLoad >= targetLoad ) ++node;
        c += cinc;
      }
      cinc *= -1;  c += cinc;
      b += binc;
    }
    binc *= -1;  b += binc;
  }

  for ( int i=0; i<patchMap->numPatches(); ++i ) {
    assignedNode[i] += unusedNodes;
  }
  sortNodesAndAssign(assignedNode);
  delete [] assignedNode; 
}

//----------------------------------------------------------------------
void WorkDistrib::mapComputes(void)
{
  PatchMap *patchMap = PatchMap::Object();
  ComputeMap *computeMap = ComputeMap::Object();
  CProxy_Node nd(CkpvAccess(BOCclass_group).node);
  Node *node = nd.ckLocalBranch();

  DebugM(3,"Mapping computes\n");

  computeMap->allocateCids();

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
      if (node->simParameters->useOptPME) {
	mapComputeHomePatches(optPmeType);
	if ( node->simParameters->pressureProfileEwaldOn )
	  mapComputeHomePatches(computeEwaldType);
      }
      else {
	mapComputeHomePatches(computePmeType);
	if ( node->simParameters->pressureProfileEwaldOn )
	  mapComputeHomePatches(computeEwaldType);
      }
    }
#else
    if (node->simParameters->useOptPME) {
      mapComputeHomePatches(optPmeType);
      if ( node->simParameters->pressureProfileEwaldOn )
	mapComputeHomePatches(computeEwaldType);
    }
    else {      
      mapComputeHomePatches(computePmeType);
      if ( node->simParameters->pressureProfileEwaldOn )
	mapComputeHomePatches(computeEwaldType);
    }
#endif
  }

  if ( node->simParameters->globalForcesOn ) {
    DebugM(2,"adding ComputeGlobal\n");
    mapComputeHomePatches(computeGlobalType);
  }

  if ( node->simParameters->extForcesOn )
    mapComputeHomePatches(computeExtType);

#ifdef NAMD_CUDA
  mapComputeNode(computeNonbondedCUDAType);
  mapComputeHomeTuples(computeExclsType);
  mapComputePatch(computeSelfExclsType);
#endif

  mapComputeNonbonded();

  // If we're doing true pair interactions, no need for bonded terms.
  // But if we're doing within-group interactions, we do need them.
  if ( !node->simParameters->pairInteractionOn || 
      node->simParameters->pairInteractionSelf) { 
    mapComputeHomeTuples(computeBondsType);
    mapComputeHomeTuples(computeAnglesType);
    mapComputeHomeTuples(computeDihedralsType);
    mapComputeHomeTuples(computeImpropersType);
    mapComputeHomeTuples(computeCrosstermsType);
    mapComputePatch(computeSelfBondsType);
    mapComputePatch(computeSelfAnglesType);
    mapComputePatch(computeSelfDihedralsType);
    mapComputePatch(computeSelfImpropersType);
    mapComputePatch(computeSelfCrosstermsType);
  }

  if ( node->simParameters->drudeOn ) {
    mapComputeHomeTuples(computeTholeType);
    mapComputePatch(computeSelfTholeType);
    mapComputeHomeTuples(computeAnisoType);
    mapComputePatch(computeSelfAnisoType);
  }

  if ( node->simParameters->eFieldOn )
    mapComputePatch(computeEFieldType);
  /* BEGIN gf */
  if ( node->simParameters->mgridforceOn )
    mapComputeHomePatches(computeGridForceType);
  /* END gf */
  if ( node->simParameters->stirOn )
    mapComputePatch(computeStirType);
  if ( node->simParameters->sphericalBCOn )
    mapComputePatch(computeSphericalBCType);
  if ( node->simParameters->cylindricalBCOn )
    mapComputePatch(computeCylindricalBCType);
  if ( node->simParameters->tclBCOn ) {
    mapComputeHomePatches(computeTclBCType);
  }
  if ( node->simParameters->constraintsOn )
    mapComputePatch(computeRestraintsType);
  if ( node->simParameters->consForceOn )
    mapComputePatch(computeConsForceType);
  if ( node->simParameters->consTorqueOn )
    mapComputePatch(computeConsTorqueType);

    // store the latest compute map
  SimParameters *simParams = Node::Object()->simParameters;
  if (simParams->storeComputeMap) {
    computeMap->saveComputeMap(simParams->computeMapFilename);
  }
    // override mapping decision
  if (simParams->loadComputeMap) {
    computeMap->loadComputeMap(simParams->computeMapFilename);
    CkPrintf("ComputeMap has been loaded from %s.\n", simParams->computeMapFilename);
  }
}

//----------------------------------------------------------------------
void WorkDistrib::mapComputeHomeTuples(ComputeType type)
{
  PatchMap *patchMap = PatchMap::Object();
  ComputeMap *computeMap = ComputeMap::Object();
  CProxy_Node nd(CkpvAccess(BOCclass_group).node);
  Node *node = nd.ckLocalBranch();

  int numNodes = node->numNodes();
  int numPatches = patchMap->numPatches();
  ComputeID cid;
  PatchIDList basepids;

  for(int i=0; i<numNodes; i++) {
    patchMap->basePatchIDList(i,basepids);
    if ( basepids.size() ) {
      cid=computeMap->storeCompute(i,basepids.size(),type);
      for(int j=0; j<basepids.size(); ++j) {
        patchMap->newCid(basepids[j],cid);
        computeMap->newPid(cid,basepids[j]);
      }
    }
  }
}

//----------------------------------------------------------------------
void WorkDistrib::mapComputeHomePatches(ComputeType type)
{
  PatchMap *patchMap = PatchMap::Object();
  ComputeMap *computeMap = ComputeMap::Object();
  CProxy_Node nd(CkpvAccess(BOCclass_group).node);
  Node *node = nd.ckLocalBranch();

  int numNodes = node->numNodes();
  int numPatches = patchMap->numPatches();
  ComputeID *cid = new ComputeID[numNodes];

  for(int i=0; i<numNodes; i++) {
    if ( patchMap->numPatchesOnNode(i) ) {
      cid[i]=computeMap->storeCompute(i,patchMap->numPatchesOnNode(i),type);
    }
  }

  PatchID j;

  for(j=0;j<numPatches;j++)
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


//----------------------------------------------------------------------
void WorkDistrib::mapComputeNode(ComputeType type)
{
  PatchMap *patchMap = PatchMap::Object();
  ComputeMap *computeMap = ComputeMap::Object();

  PatchID i;
  ComputeID cid;

  for(int i=0; i<CkNumPes(); i++) {
    computeMap->storeCompute(i,0,type);
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
  CProxy_Node nd(CkpvAccess(BOCclass_group).node);
  Node *node = nd.ckLocalBranch();

  PatchID oneAway[PatchMap::MaxOneOrTwoAway];
  PatchID oneAwayTrans[PatchMap::MaxOneOrTwoAway];

  PatchID i;
  ComputeID cid;
  int numNeighbors;
  int j;

  for(i=0; i<patchMap->numPatches(); i++) // do the self 
  {
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////^M
//// Osman Sarood^M
//// Parallel Input Change^M
//// Comment some lines^M
 int64 numAtoms = patchMap->patch(i)->getNumAtoms();
#ifdef  MEM_OPT_VERSION
    int64 numFixed = 0;  // avoid overflow^M
#else
        int64 numFixed = patchMap->patch(i)->getNumFixedAtoms();  // avoid overflow^M
#endif
//////////////////////////////////////////////////
    int numPartitions = 0;
    int divide = node->simParameters->numAtomsSelf;
    if (divide > 0) {
      numPartitions = (int) ( 0.5 +
        (numAtoms*numAtoms-numFixed*numFixed) / (double)(2*divide*divide) );
    }
    if (numPartitions < 1) numPartitions = 1;
    if ( numPartitions > node->simParameters->maxSelfPart )
			numPartitions = node->simParameters->maxSelfPart;
    // self-interaction
    DebugM(4,"Mapping " << numPartitions << " ComputeNonbondedSelf objects for patch " << i << "\n");
//    iout <<"Self numPartitions = " <<numPartitions <<" numAtoms " <<numAtoms <<std::endl;
    for(int partition=0; partition < numPartitions; partition++)
    {
      cid=computeMap->storeCompute(patchMap->node(i),1,
				   computeNonbondedSelfType,
				   partition,numPartitions);
      computeMap->newPid(cid,i);
      patchMap->newCid(i,cid);
    }
  }

  for(int p1=0; p1 <patchMap->numPatches(); p1++) // do the pairs
  {
    // this only returns half of neighbors, which is what we want
    numNeighbors=patchMap->oneOrTwoAwayNeighbors(p1,oneAway,oneAwayTrans);
    for(j=0;j<numNeighbors;j++)
    {
	int p2 = oneAway[j];
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////^M
//// Osman Sarood^M
//// Parallel Input Change^M
//// Comment some lines^M
        int64 numAtoms1 = patchMap->patch(p1)->getNumAtoms();
        int64 numAtoms2 = patchMap->patch(p2)->getNumAtoms();

#ifdef  MEM_OPT_VERSION
        int64 numFixed1 = patchMap->patch(p1)->getNumFixedAtoms();
        int64 numFixed2 = patchMap->patch(p2)->getNumFixedAtoms();
#else

        int64 numFixed1 = 0;
        int64 numFixed2 = 0;
#endif

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

        const int nax = patchMap->numaway_a();  // 1 or 2
        const int nay = patchMap->numaway_b();  // 1 or 2
        const int naz = patchMap->numaway_c();  // 1 or 2
        const int ia1 = patchMap->index_a(p1);
        const int ia2 = patchMap->index_a(p2);
        const int ib1 = patchMap->index_b(p1);
        const int ib2 = patchMap->index_b(p2);
        const int ic1 = patchMap->index_c(p1);
        const int ic2 = patchMap->index_c(p2);
	int distance = 3;
 	if ( ia1 == ia2 ) --distance;
 	else if ( ia1 == ia2 + nax - 1 ) --distance;
 	else if ( ia1 + nax - 1 == ia2 ) --distance;
 	if ( ib1 == ib2 ) --distance;
 	else if ( ib1 == ib2 + nay - 1 ) --distance;
 	else if ( ib1 + nay - 1 == ib2 ) --distance;
 	if ( ic1 == ic2 ) --distance;
 	else if ( ic1 == ic2 + naz - 1 ) --distance;
 	else if ( ic1 + naz - 1 == ic2 ) --distance;
	int divide = 0;
	if ( distance == 0 ) {
	  divide = node->simParameters->numAtomsSelf2;
        } else if (distance == 1) {
	  divide = node->simParameters->numAtomsPair;
	} else {
	  divide = node->simParameters->numAtomsPair2;
	}
        int numPartitions = 0;
	if (divide > 0) {
          numPartitions = (int) ( 0.5 +
	    (numAtoms1*numAtoms2-numFixed1*numFixed2)/(double)(divide*divide) );
	}
        if ( numPartitions < 1 ) numPartitions = 1;
        if ( numPartitions > node->simParameters->maxPairPart )
			numPartitions = node->simParameters->maxPairPart;
//	if ( numPartitions > 1 ) iout << "Mapping " << numPartitions << " ComputeNonbondedPair objects for patches " << p1 << "(" << numAtoms1 << ") and " << p2 << "(" << numAtoms2 << ")\n" << endi;
	for(int partition=0; partition < numPartitions; partition++)
	{
	  cid=computeMap->storeCompute(
		patchMap->basenode(patchMap->downstream2(p1,p2)),
		2,computeNonbondedPairType,partition,numPartitions);
	  computeMap->newPid(cid,p1);
	  computeMap->newPid(cid,p2,oneAwayTrans[j]);
	  patchMap->newCid(p1,cid);
	  patchMap->newCid(p2,cid);
        }
    }
  }

}

//----------------------------------------------------------------------
void WorkDistrib::messageEnqueueWork(Compute *compute) {
  LocalWorkMsg *msg = compute->localWorkMsg;
  int seq = compute->sequence();

  if ( seq < 0 ) {
    NAMD_bug("compute->sequence() < 0 in WorkDistrib::messageEnqueueWork");
  } else {
    SET_PRIORITY(msg,seq,compute->priority());
  }

  msg->compute = compute; // pointer is valid since send is to local Pe
  int type = compute->type();

  CProxy_WorkDistrib wdProxy(CkpvAccess(BOCclass_group).workDistrib);
  switch ( type ) {
  case computeExclsType:
  case computeSelfExclsType:
    wdProxy[CkMyPe()].enqueueExcls(msg);
    break;
  case computeBondsType:
  case computeSelfBondsType:
    wdProxy[CkMyPe()].enqueueBonds(msg);
    break;
  case computeAnglesType:
  case computeSelfAnglesType:
    wdProxy[CkMyPe()].enqueueAngles(msg);
    break;
  case computeDihedralsType:
  case computeSelfDihedralsType:
    wdProxy[CkMyPe()].enqueueDihedrals(msg);
    break;
  case computeImpropersType:
  case computeSelfImpropersType:
    wdProxy[CkMyPe()].enqueueImpropers(msg);
    break;
  case computeTholeType:
  case computeSelfTholeType:
    wdProxy[CkMyPe()].enqueueThole(msg);
    break;
  case computeAnisoType:
  case computeSelfAnisoType:
    wdProxy[CkMyPe()].enqueueAniso(msg);
    break;
  case computeCrosstermsType:
  case computeSelfCrosstermsType:
    wdProxy[CkMyPe()].enqueueCrossterms(msg);
    break;
  case computeNonbondedSelfType:
    switch ( seq % 2 ) {
    case 0:
      wdProxy[CkMyPe()].enqueueSelfA(msg);
      break;
    case 1:
      wdProxy[CkMyPe()].enqueueSelfB(msg);
      break;
    default:
      NAMD_bug("WorkDistrib::messageEnqueueSelf case statement error!");
    }
    break;
  case computeNonbondedPairType:
    switch ( seq % 2 ) {
    case 0:
      wdProxy[CkMyPe()].enqueueWorkA(msg);
      break;
    case 1:
      wdProxy[CkMyPe()].enqueueWorkB(msg);
      break;
    case 2:
      wdProxy[CkMyPe()].enqueueWorkC(msg);
      break;
    default:
      NAMD_bug("WorkDistrib::messageEnqueueWork case statement error!");
    }
    break;
  case computeNonbondedCUDAType:
#ifdef NAMD_CUDA
    // CkPrintf("CUDA %d %d %x\n", CkMyPe(), seq, compute->priority());
    wdProxy[CkMyPe()].enqueueCUDA(msg);
#else
    msg->compute->doWork();
#endif
    break;
  case computePmeType:
    // CkPrintf("PME %d %d %x\n", CkMyPe(), seq, compute->priority());
#ifdef NAMD_CUDA
    wdProxy[CkMyPe()].enqueuePme(msg);
#else
    msg->compute->doWork();
#endif
    break;
  case optPmeType:
    // CkPrintf("PME %d %d %x\n", CkMyPe(), seq, compute->priority());
#ifdef NAMD_CUDA
    wdProxy[CkMyPe()].enqueuePme(msg);
#else
    msg->compute->doWork();
#endif
    break;
  default:
    wdProxy[CkMyPe()].enqueueWork(msg);
  }
}

void WorkDistrib::enqueueWork(LocalWorkMsg *msg) {
  msg->compute->doWork();
  if ( msg->compute->localWorkMsg != msg )
    NAMD_bug("WorkDistrib LocalWorkMsg recycling failed!");
}

void WorkDistrib::enqueueExcls(LocalWorkMsg *msg) {
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

void WorkDistrib::enqueueThole(LocalWorkMsg *msg) {
  msg->compute->doWork();
  if ( msg->compute->localWorkMsg != msg )
    NAMD_bug("WorkDistrib LocalWorkMsg recycling failed!");
}

void WorkDistrib::enqueueAniso(LocalWorkMsg *msg) {
  msg->compute->doWork();
  if ( msg->compute->localWorkMsg != msg )
    NAMD_bug("WorkDistrib LocalWorkMsg recycling failed!");
}

void WorkDistrib::enqueueCrossterms(LocalWorkMsg *msg) {
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

void WorkDistrib::enqueueCUDA(LocalWorkMsg *msg) {
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
  SimParameters *simParams = Node::Object()->simParameters;
  Bool lesOn = simParams->lesOn;
  Random vel_random(simParams->randomSeed);

  int lesReduceTemp = lesOn && simParams->lesReduceTemp;
  BigReal tempFactor = lesReduceTemp ? 1.0 / simParams->lesFactor : 1.0;

  kbT = Temp*BOLTZMAN;

  //  Loop through all the atoms and assign velocities in
  //  the x, y and z directions for each one
  for (i=0; i<totalAtoms; i++)
  {
    if (structure->atommass(i) <= 0.) {
      kbToverM = 0.;
    } else {
      kbToverM = sqrt(kbT *
        ( lesOn && structure->get_fep_type(i) ? tempFactor : 1.0 ) /
			  structure->atommass(i) );
    }

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

#if USE_TOPOMAP 

//Specifically designed for BGL and other 3d Tori architectures
//Partition Torus and Patch grid together using recursive bisection.
int WorkDistrib::assignPatchesTopoGridRecBisection() {
  
  PatchMap *patchMap = PatchMap::Object();
  int *assignedNode = new int[patchMap->numPatches()];
  int numNodes = Node::Object()->numNodes();
  SimParameters *simParams = Node::Object()->simParameters;
  int usedNodes = numNodes;
  
  if ( simParams->noPatchesOnZero && numNodes > 1 ) usedNodes -= 1;
  RecBisection recBisec(patchMap->numPatches(), PatchMap::Object());
  
  int xsize = 0, ysize = 0, zsize = 0;
  
  // Right now assumes a T*** (e.g. TXYZ) mapping
  TopoManager tmgr;
  xsize = tmgr.getDimNX();
  ysize = tmgr.getDimNY();
  zsize = tmgr.getDimNZ();
  
  //Fix to not assign patches to processor 0
  int rc = recBisec.partitionProcGrid(xsize, ysize, zsize, assignedNode);
 
  delete [] assignedNode;

  return rc;
}
#endif

#include "WorkDistrib.def.h"

