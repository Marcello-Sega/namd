#include <stdio.h>
#include "BOCgroup.h"
#include "ParallelIOMgr.decl.h"
#include "ParallelIOMgr.h"
#include "Molecule.h"
#include "Node.h"
#include "Node.decl.h"
#include "NamdState.h"
#include "WorkDistrib.h"
#include "PDB.h"
#include "PatchMap.h"
#include "PatchMap.inl"
#include "packmsg.h"
#include "HomePatch.h"
#ifdef MEM_OPT_VERSION

// Method used to distribute assigned nodes inof
void ParallelIOMgr::initAssignedNodesParallel()
{
  PatchMap *patchMap=Node::Object()->getPatchMap();
  assignedNodesParallel=new int[patchMap->numPatches()];
  for(int p=0;p<patchMap->numPatches();p++)
  {
          assignedNodesParallel[p]=patchMap->node(p);
  }
  broadcastPatchInfo();
}

// Helper method which calls method to calculate atom-to-patch assignment method
void ParallelIOMgr::calcAtomsEachPatch()
{
  CProxy_ParallelIOMgr pIO(thisgroup);
  pIO.caclNumAtomsInEachPatchParallel();
}
#endif
// does the atom-to-patch assignment
void ParallelIOMgr::caclNumAtomsInEachPatchParallel(){
#ifdef MEM_OPT_VERSION
  PatchMap *patchMap = PatchMap::Object();
  int numPatches=patchMap->numPatches();
  if(Node::Object()->ioMgr->isInputProc(CkMyPe()))
  {
    patchMap->initTmpPatchAtomsList();
    StringList *current;
    int i;
    CProxy_Node nd(CpvAccess(BOCclass_group).node);
    Node *node = nd.ckLocalBranch();
    CProxy_PatchMgr pm(CpvAccess(BOCclass_group).patchMgr);
    PatchMgr *patchMgr = pm.ckLocalBranch();
    SimParameters *params = node->simParameters;
    Molecule *molecule = node->molecule;
    PDB *pdb = node->pdb;


    int numAtomsPerProc=Node::Object()->ioMgr->getNewTotalAtoms();

    vector<int> *eachPatchAtomList = patchMap->getTmpPatchAtomsList();
    int *patchAtomCnt = new int[numPatches];    
    memset(patchAtomCnt, 0, sizeof(int)*numPatches);

    const Lattice lattice = params->lattice;

    Position eachAtomPos;
    if (params->splitPatch == SPLIT_PATCH_HYDROGEN)
    {

      int aid, pid=0;
      for(i=0; i < numAtomsPerProc; i++)
        {        
        aid = hydrogenGroupPar[i].atomID;
        Node::Object()->ioMgr->get_position_for_atom_Parallel(&eachAtomPos,aid);
        if(hydrogenGroupPar[i].isMP)
            pid = patchMap->assignToPatch(eachAtomPos,lattice);
        
	patchAtomCnt[pid]++;
        eachPatchAtomList[pid].push_back(aid);
        }
      }
    else
    {
      for(i=0; i < numAtomsPerProc; i++)
      {
	Node::Object()->ioMgr->get_position_for_atom_Parallel(&eachAtomPos,i);
        int pid ;
        pid= patchMap->assignToPatch(eachAtomPos,lattice);        
        patchAtomCnt[pid]++;
        eachPatchAtomList[pid].push_back(i);
        }
    }   
      AtomsPerPatchMsg *msg1=new AtomsPerPatchMsg;
      msg1->patchAtomCnt=patchAtomCnt;
      msg1->numPatches=numPatches;

      CProxy_ParallelIOMgr pIO(thisgroup);
      pIO[0].receiveNumAtomsPerPatch((AtomsPerPatchMsg*)AtomsPerPatchMsg::pack(msg1));
  }
#endif
}

#ifdef MEM_OPT_VERSION
//gets the velocities for a given atom id
void ParallelIOMgr::get_position_for_patch_creationVel(Vector *pos,int aid,int aidIdx)
{
  int i=aidIdx;
  FullAtom *recvAt=this->finalAtomList.begin();
  pos->x=recvAt[i].velocity.x;
  pos->y=recvAt[i].velocity.y;
  pos->z=recvAt[i].velocity.z;
}

// gets the position for a given atoms id
void ParallelIOMgr::get_position_for_patch_creation(Vector *pos,int aid,int aidIdx)
{
  int i=aidIdx;
  FullAtom *recvAt=this->finalAtomList.begin();
  pos->x=recvAt[i].position.x;
  pos->y=recvAt[i].position.y;
  pos->z=recvAt[i].position.z;
}

// helper method which broadcasts the patch-to-node assignment
void ParallelIOMgr::broadcastPatchInfo()
{
  CProxy_ParallelIOMgr pIO(thisgroup);
  PatchMap *patchMap=Node::Object()->getPatchMap();
  for(int p=0;p<CkNumPes();p++)
  {
    NodeAssignmentMsg *msg=new NodeAssignmentMsg;
    msg->totalPatches=patchMap->numPatches();
    msg->atomsPerPatchParallel=atomsPerPatchParallel;
    msg->nodeAssigned=assignedNodesParallel;                
    pIO[p].sendAssignedNode((NodeAssignmentMsg*)NodeAssignmentMsg::pack(msg));
  }  
}
#endif

// entry method through which each chare obtains patch-to-node assignment
void ParallelIOMgr::sendAssignedNode(NodeAssignmentMsg *msg)
{
#ifdef MEM_OPT_VERSION
  CProxy_Node nd(CpvAccess(BOCclass_group).node);
  Node *node = nd.ckLocalBranch();
  Molecule *molecule = node->molecule;

  PatchMap *patchMap = PatchMap::Object();


  if(patchMap->numPatches()==0) NAMD_bug("Number of patches are zero in WorkDistrib.C!");
    assignedNodesParallel=new int[patchMap->numPatches()];
  for(int i=0;i<patchMap->numPatches();i++)
  {
    assignedNodesParallel[i]=msg->nodeAssigned[i];
  }
  atomsPerPatchParallel=msg->atomsPerPatchParallel;
  Node::Object()->ioMgr->createAtomHoldPatchTransfer(msg->totalPatches);
#endif
}

#ifdef MEM_OPT_VERSION
// to free space taken by patches
void ParallelIOMgr::deleteAtomHoldPatchTransfer()
{
  delete [] patchesAfterPatchExchange;
  delete [] patchAtomList;// = new vector<int>[numPatches];
  delete [] patchAtomIdMappingList;// = new vector<int>[numPatches];
  hydrogenGroupPerProc.resize(0);
}

// space allocated to hold atom data just before the creation of patches
void ParallelIOMgr::createAtomHoldPatchTransfer(int numPatches)
{

  currentAtomsAfterPatchExchange=0;
  int totalAtomsOnProc=0;
  for(int i=0;i<numPatches;i++)
  {
                
    if(CkMyPe()==assignedNodesParallel[i])
    {
               
      totalAtomsOnProc+=atomsPerPatchParallel[i];
    }       
  }
  patchesAfterPatchExchange=new int[totalAtomsOnProc];
  patchAtomList = new vector<int>[numPatches];
  patchAtomIdMappingList = new vector<int>[numPatches];
  hydrogenGroupPerProc.resize(totalAtomsOnProc);
}

// method that creates the patches after the data has moved to approporate processors
void ParallelIOMgr::createPatchesThisProc()
{
  Vector *velocities = new Velocity[currentAtomsAfterPatchExchange];
  if ( simParameters->initialTemp < 0.0 ) {
  }
  else
    Node::Object()->workDistrib->random_velocities_parallel(simParameters->initialTemp,molecule,velocities,currentAtomsAfterPatchExchange);
  hydrogenGroupPerProc.sort();
  mappingHydGPPerProc=new int[hydrogenGroupPerProc.size()];
  mapHydGPPerProc();
  for(int i=0;i<hydrogenGroupPerProc.size();i++)
  {
    int aid=hydrogenGroupPerProc[i].atomID;
    int aidIdx=hydrogenGroupPerProc[i].atomHoldIdxPar;
    int patchID=patchesAfterPatchExchange[aidIdx];
    patchAtomList[patchID].push_back(aid);
    patchAtomIdMappingList[patchID].push_back(aidIdx);
  }
  for(int i=0; i < Node::Object()->getPatchMap()->numPatches(); i++)
  {
    if(assignedNodesParallel[i]==CkMyPe())
    {
      FullAtomList *onePatchAtoms = new FullAtomList;
      Node::Object()->workDistrib->fillOnePatchCreationParallelIO(i, onePatchAtoms, velocities);
      Node::Object()->getPatchMgr()->createHomePatch(i,*onePatchAtoms);
    }
  }
  delete [] velocities;
}
#endif
// entry method that receives the atoms transferred in order to complete migration groups at each proc
void ParallelIOMgr::receiveMigrationChildren(MoveFullAtomsMsg *msg)
{
#ifdef MEM_OPT_VERSION

  atomsEnteringProc=msg->atom.size();
  FullAtom *recvList=msg->atom.begin();
  for(int i=0;i<msg->atom.size();i++)
  {

    FullAtom newAt;
    newAt=recvList[i];
    recvAtomList.add(newAt);
  }

  if(checkForSendAndRecv() && hydrogenGPUpated==false)
    hydrogenGPUpated=true;
    
#endif
}

// used to recevie the number of atoms per patch so that the total can be calculated. Called on pe0
void ParallelIOMgr::receiveNumAtomsPerPatch(AtomsPerPatchMsg *msg1)
{
#ifdef MEM_OPT_VERSION
  procsReceived++;
  numProcsSentNumAtomsPerPatch++;
  for(int i=0;i<msg1->numPatches;i++)
    atomsPerPatchParallel[i]=atomsPerPatchParallel[i]+msg1->patchAtomCnt[i];
  
  PatchMap *patchMap=Node::Object()->getPatchMap();
  int numPatches = patchMap->numPatches();
  int numAts=0;
  for(int i=0; i < numPatches; i++)
    numAts+=atomsPerPatchParallel[i];
    
  if(procsReceived==numInputProcs)
  {
    for(int i=0;i<patchMap->numPatches();i++)
    {
      HomePatch *patch = new HomePatch(i, atomsPerPatchParallel[i]);
      patchMap->registerMyPatch(i,patch);
    }
  }
  delete [] msg1->patchAtomCnt;
  delete msg1;
#endif
}

// helper function to intiate final atom transfer before patch creation
void ParallelIOMgr::atomExchangeForPatchCreation()
{
#ifdef MEM_OPT_VERSION
    CProxy_ParallelIOMgr pIO(thisgroup);
    for(int i=0;i<numInputProcs;i++)
      pIO[inputProcArray[i]].sendCoordinatesToProcs();  
#endif
}

// preparing messages to be sent to other procs so that patches can be created
void ParallelIOMgr::sendCoordinatesToProcs()
{
#ifdef MEM_OPT_VERSION

  AtomListForPatchCreateMsg **msg=new AtomListForPatchCreateMsg*[CkNumPes()];
  FullAtomList *migrationAtoms=new FullAtomList[CkNumPes()];


  PatchMap *patchMap = Node::Object()->getPatchMap();        
  if(isInputProc(CkMyPe()))
  {
    int **atomsPerPatch=new int*[CkNumPes()];
    int **patcheNumsToSend=new int*[CkNumPes()];
    int *patchesToSend=new int[CkNumPes()];

    for(int inProc=0;inProc<CkNumPes();inProc++)
    {
      int patchesToThisProc=0;
      for(int i=0;i<patchMap->numPatches();i++)
      {
	if(inProc==assignedNodesParallel[i])
	  patchesToThisProc++;
	             
      }
      atomsPerPatch[inProc]=new int[patchesToThisProc];
      patcheNumsToSend[inProc]=new int[patchesToThisProc];
      patchesToSend[inProc]=0;
    }

    fflush(stdout);
    for(int p=0;p<patchMap->numPatches();p++)
    {
      vector<int> *eachPatchAtomsList = patchMap->getTmpPatchAtomsList();
      vector<int> *thisPatchAtomsList = &eachPatchAtomsList[p];
      int inProc=assignedNodesParallel[p];
      {
        HydrogenGroupID *hg = molecule->hydrogenGroup.begin();
        HydrogenGroupID *hgPar = hydrogenGroupPar.begin();

        int startAtom=getRecordOffset();

        int endAtom=startAtom+getAtomsAssignedThisProc();
        for(int at=0;at<thisPatchAtomsList->size();at++)
        {
	  int aid = thisPatchAtomsList->at(at);
	  Position pos;
	  Position vel;

	  get_position_for_atom_Parallel(&pos,aid);
	  if ( simParameters->initialTemp < 0.0 ) {
	    get_position_for_atom_ParallelVel(&vel, aid);
	  }
	  FullAtom newAt;
	  newAt.id=aid;
	  newAt.position.x=pos.x;
	  newAt.position.y=pos.y;
	  newAt.position.z=pos.z;
          if ( simParameters->initialTemp < 0.0 ) {
	    newAt.velocity.x=vel.x;
	    newAt.velocity.y=vel.y;
	    newAt.velocity.z=vel.z;
          }
          bool found=false;
          if(!isAtomSent(aid))
          {
            if(aid>=startAtom && aid<endAtom)
            {
              found=true;
	      newAt.GPID=hg[aid-recordOffset].GPID;
	      newAt.sortVal=hg[aid-recordOffset].sortVal;
	      newAt.atomsInGroup=hg[aid-recordOffset].atomsInGroup;
	      newAt.isMP=hg[aid-recordOffset].isMP;
	      newAt.migrationGroupSize=hg[aid-recordOffset].atomsInMigrationGroup;
	      newAt.MPID=hg[aid-recordOffset].MPID;
	      newAt.waterVal=hg[aid-recordOffset].waterVal;

            }
          }
          if(newRecvAtom(aid))
          {
	    for(int i=0;i<getNewTotalAtoms();i++)
            {
	      if(hgPar[i].atomID==aid)
              {
		found=true;
		newAt.GPID=hgPar[i].GPID;
		newAt.sortVal=hgPar[i].sortVal;
		newAt.atomsInGroup=hgPar[i].atomsInGroup;
		newAt.isMP=hgPar[i].isMP;
		newAt.migrationGroupSize=hgPar[i].atomsInMigrationGroup;
		newAt.MPID=hgPar[i].MPID;
		newAt.waterVal=hgPar[i].waterVal;
              }
	    }
          }
          if(found==false)
            NAMD_bug("In workDistrib.C there is a problem in sendAssignedNodeToEachInputProc()");
          #ifdef MEM_OPT_VERSION
	  newAt.sigId=atomsigIDPar(aid);
          #endif
	  newAt.exclId=atomexclsigIDPar(aid);
	  newAt.vdwType=atomvdwtypePar(aid);
	  newAt.mass=atommassPar(aid);
	  newAt.charge=atomchargePar(aid);
                                                
	  migrationAtoms[inProc].add(newAt);                                                
        }
        atomsPerPatch[inProc][patchesToSend[inProc]]=thisPatchAtomsList->size();
        patcheNumsToSend[inProc][patchesToSend[inProc]]=p;
        patchesToSend[inProc]++;


                                
      }
    }
    CProxy_Node cm(thisgroup);
    for(int inProc=0;inProc<CkNumPes();inProc++)
    {
      if(migrationAtoms[inProc].size()>0)
      {
	msg[inProc] = new AtomListForPatchCreateMsg();
	msg[inProc]->patchesToSend=patchesToSend[inProc];
	msg[inProc]->atom=migrationAtoms[inProc];
	msg[inProc]->atomsPerPatch=atomsPerPatch[inProc];
	msg[inProc]->patcheNumsToSend=patcheNumsToSend[inProc];
                        
        CProxy_ParallelIOMgr ioMgrP(thisgroup);
	ioMgrP[inProc].getPatchCoordinates((AtomListForPatchCreateMsg*)AtomListForPatchCreateMsg::pack(msg[inProc]));

      }
      delete [] atomsPerPatch[inProc];
      delete [] patcheNumsToSend[inProc];
    }
    delete [] atomsPerPatch;
    delete [] patcheNumsToSend;


  }
                
  fflush(stdout);
#endif
}
                        

// method which is called to transfer atom data before patch creation
void ParallelIOMgr::getPatchCoordinates(AtomListForPatchCreateMsg *msg)
{
#ifdef MEM_OPT_VERSION
  int totalAtomsOnProc=0;
  PatchMap *patchMap=Node::Object()->getPatchMap();
  for(int i=0;i<patchMap->numPatches();i++)
  {
    // if the patch is assigned to this proc^M
    if(CkMyPe()==assignedNodesParallel[i])
    {                        //add the num of atoms in that patch to get the total number of atoms on this proc.
     totalAtomsOnProc+=atomsPerPatchParallel[i];
    }
  }
  // for every patch in the system
  int allAtomCount=0;
  // iterate for numPatchesSent
  for(int p=0;p<msg->patchesToSend;p++)
  {
    int atomsInThisPatch=msg->atomsPerPatch[p];
    // for each of the atoms present in the patch.
    FullAtom *recvAtom=msg->atom.begin();
    for(int at=0;at<atomsInThisPatch;at++)
    {
      FullAtom newAt;
      int aid = recvAtom[allAtomCount].id;
      if ( simParameters->initialTemp < 0.0 ) {

	newAt.velocity.x=recvAtom[allAtomCount].velocity.x;
	newAt.velocity.y=recvAtom[allAtomCount].velocity.y;
	newAt.velocity.z=recvAtom[allAtomCount].velocity.z;
      }

      newAt.position.x=recvAtom[allAtomCount].position.x;
      newAt.position.y=recvAtom[allAtomCount].position.y;
      newAt.position.z=recvAtom[allAtomCount].position.z;

      patchesAfterPatchExchange[currentAtomsAfterPatchExchange]=msg->patcheNumsToSend[p];
      hydrogenGroupPerProc[currentAtomsAfterPatchExchange].atomID=aid;
      hydrogenGroupPerProc[currentAtomsAfterPatchExchange].atomsInGroup=recvAtom[allAtomCount].atomsInGroup;
      hydrogenGroupPerProc[currentAtomsAfterPatchExchange].GPID=recvAtom[allAtomCount].GPID;
      hydrogenGroupPerProc[currentAtomsAfterPatchExchange].sortVal=recvAtom[allAtomCount].sortVal;
      hydrogenGroupPerProc[currentAtomsAfterPatchExchange].isGP=1;
      hydrogenGroupPerProc[currentAtomsAfterPatchExchange].isMP=recvAtom[allAtomCount].isMP;
      hydrogenGroupPerProc[currentAtomsAfterPatchExchange].atomsInMigrationGroup=recvAtom[allAtomCount].migrationGroupSize;
      hydrogenGroupPerProc[currentAtomsAfterPatchExchange].MPID=recvAtom[allAtomCount].MPID;
      hydrogenGroupPerProc[currentAtomsAfterPatchExchange].waterVal=recvAtom[allAtomCount].waterVal;

      if(hydrogenGroupPerProc[currentAtomsAfterPatchExchange].atomsInGroup==0) 
      hydrogenGroupPerProc[currentAtomsAfterPatchExchange].isGP = 0;
      hydrogenGroupPerProc[currentAtomsAfterPatchExchange].atomHoldIdxPar=currentAtomsAfterPatchExchange;

      newAt.id=aid;
      newAt.mass=recvAtom[allAtomCount].mass;
      newAt.charge=recvAtom[allAtomCount].charge;
      newAt.sigId=recvAtom[allAtomCount].sigId;
      newAt.exclId=recvAtom[allAtomCount].exclId;
      newAt.vdwType=recvAtom[allAtomCount].vdwType;
      finalAtomList.add(newAt);

      currentAtomsAfterPatchExchange++;
      allAtomCount++;
    }
  }
  delete [] msg->atomsPerPatch;
  delete [] msg->patcheNumsToSend;
  delete msg;
#endif
}
                
#ifdef MEM_OPT_VERSION
// initialization to determine which procs would be the input procs
void ParallelIOMgr::initInputProcArray(int inputProcs)
{
  inputProcArray=new int[inputProcs];
  numInputProcs=inputProcs;
  for(int p=0;p<inputProcs;p++)
    inputProcArray[p]=p;
}

// to determine whether a given proc is an input proc
bool ParallelIOMgr::isInputProc(int proc)
{
  for(int i=0;i<numInputProcs;i++)
  {
    if(proc==inputProcArray[i])
      return true;
  }
  return false;
}

//get the idx for a proc in the input proc array
int ParallelIOMgr::inputProcIndex(int proc)
{
  for(int i=0;i<numInputProcs;i++)
  {
    if(proc==inputProcArray[i])
      return i;
  }
  return -1;
}

// checking for conditions not yet compatible with parallel IO
void ParallelIOMgr::checkConfigParas()
{
  if(Node::Object()->simParameters->numinputprocs>CkNumPes())
    NAMD_bug("Num of Input Processors exceed the total number of processors! They can not be more than p-1.\n");
  
  if(Node::Object()->simParameters->langevinOn)
    NAMD_bug("Langevin On can not work with ParallelIO\n");

  if(Node::Object()->simParameters->pairInteractionOn)
    NAMD_bug("pairInteractionOn can not work with ParallelIO\n");

  if(Node::Object()->simParameters->alchFepOn)
    NAMD_bug("alchFepOn can not work with ParallelIO\n");

  if(Node::Object()->simParameters->alchThermIntOn)
    NAMD_bug("alchThermIntOn can not work with ParallelIO\n");

  if(Node::Object()->simParameters->lesOn)
    NAMD_bug("lesOn can not work with ParallelIO\n");

  if((Node::Object()->state->getConfigList()->find("numinputprocs"))==NULL)
    NAMD_bug("Please specify numinputprocs in config file for Parallel IO to use\n");

  if((Node::Object()->state->getConfigList()->find("bincoordinates"))==NULL)
    NAMD_bug("Please specify the bincoordinates file for ParallelIO\n");
  
}

// initialization needed for the parallel IO manager. Entry point through Scripttcl
void ParallelIOMgr::initParallelInput(){
  checkConfigParas();
  CProxy_Node cml(thisgroup);
  setPointers(Node::Object()->molecule,cml,Node::Object()->state,Node::Object()->getPatchMgr(),Node::Object()->simParameters);
  int inputProcs=-1;
  inputProcs=atoi((state->getConfigList()->find("numinputprocs"))->data);

for(int i=0;i<CkNumPes();i++)
{
  ConfigListMessage *pmsg=new ConfigListMessage();
  if(state->getConfigList()!=NULL) pmsg->cfgList=state->getConfigList();
  initInputProcArray(inputProcs);
  pmsg->numAtoms=molecule->numAtoms;
  pmsg->doFlip=doFlip;
  CProxy_ParallelIOMgr pIO(thisgroup);
  pIO[i].readCoordinatesAndVel((ConfigListMessage*)ConfigListMessage::pack(pmsg));
}

}
#endif
// read assigned atoms positions and velocities
void ParallelIOMgr::readCoordinatesAndVel(ConfigListMessage *pmsg)
{
#ifdef MEM_OPT_VERSION
  int inP=-1;
  inP=atoi((pmsg->cfgList->find("numinputprocs"))->data);

  if(CkMyPe())
  {
    state=new NamdState();
    Node::Object()->state=state;
    state->setConfigList(pmsg->cfgList);
  }
  initInputProcArray(inP);
  totalAtomsSys=pmsg->numAtoms;
  if(isInputProc(CkMyPe()))
  {
    char *coorFileN;

    if(state==NULL) printf("PROC#%d There is a PROBLEM IN NODE.C STATE is NULL !!!!!!!!!!!!!!!!!!!!!!!!!!!!\n\n",CkMyPe());
    if (state->getConfigList()==NULL) printf("PROC#%d There is a PROBLEM IN NODE.C CONFIGLIST is NULL !!!!!!!!!!!!!!!!!!!!!!!!!!!!\n\n",CkMyPe());
    coorFileN=(state->getConfigList()->find("bincoordinates"))->data;
    int inputProcIdx=inputProcIndex(CkMyPe());
    int numInProcs=inP;
    int numAtom=pmsg->numAtoms;
    bool flip=pmsg->doFlip;
    initParallelIOMgr(CkMyPe(),numInProcs,numAtom,coorFileN,flip,inputProcArray);
    numInputProcs=numInProcs;
    int atomsAssignedThisProc=getAtomsAssignedThisProc();
    if(CkMyPe()==0)
    {
      this->setAtomsAssignedThisProc(atomsAssignedThisProc);
      this->setRecordOffset(getAtomsPerProcOnetoNMinus1()*inputProcIdx);
    }
    double startT=CmiWallTimer();

    readCoordinates();
    StringList *cfgFile=(state->getConfigList()->find("parameters"));
    char *fileN;
    if(state->getConfigList()->find("binvelocities")==NULL)
      fileN=NULL;
    else
    {
      fileN=(state->getConfigList()->find("binvelocities"))->data;
      this->readAtomInfoParallelInputVel(fileN,CkMyPe());
    }
  }

  delete pmsg;
#endif
}
#ifdef MEM_OPT_VERSION
// check whether an atom is recevived as a result of migration exchange
bool ParallelIOMgr::newRecvAtom(int atomID)
{
  FullAtom *recvAtom=this->recvAtomList.begin();
  for(int i=0;i<atomsEnteringProc;i++)
  {
    if(recvAtom[i].id==atomID)
      return true;
                
  }
  return false;
}

// index for the atom recvd due to migration atoms exchange
long long ParallelIOMgr::newRecvAtomIdx(int atomID)
{
  FullAtom *recvAtom=this->recvAtomList.begin();
  for(long long i=0;i<atomsEnteringProc;i++)
  {
    if(recvAtom[i].id==atomID)
      return i;
                
  }
  return -1;
}

// whether an  atom is no longer part of this proc due to migration atom exchange
bool ParallelIOMgr::isAtomSent(int atomID)
{
  for(int i=0;i<atomsLeavingProc;i++)
  {
    if(atomsLeaving[i]==atomID)
      return true;
      
  }
  return false;
}

// whether all send and recveives have oiccured for migration atoms exchange
bool ParallelIOMgr::checkForSendAndRecv()
{
  bool hasToRecv=false;
  bool hasToSend=false;
  FullAtom *migList=migrationList->begin();
  for(int i=0;i<getNumMigrationAtoms();i++)
  {
    if(migList[i].readProc==CkMyPe())
      hasToSend=true;
    if(migList[i].destProc==CkMyPe())
      hasToRecv=true;
      
  }
  if(hasToRecv && hasToSend)
  {
    if(atomsLeavingProc>0 && atomsEnteringProc>0)
      return true;
      
  }
  else if(hasToRecv && !hasToSend)
  {
    if(atomsEnteringProc>0)
      return true;
      
  }
  else if(!hasToRecv && hasToSend)
  {
    if(atomsLeavingProc>0)
      return true;
      
  }
  else
  {
    return true;
  }
  return false;
}

// build a mapping for hydrogen group
void ParallelIOMgr::mapHydGPPerProc()
{
   HydrogenGroupID *hg = hydrogenGroupPerProc.begin();
   for(int i=0;i<hydrogenGroupPerProc.size();i++)
   {
    mappingHydGPPerProc[hg[i].atomHoldIdxPar]=i;
  }
}

// value of vdwtype for a given atom id
short ParallelIOMgr::atomvdwtypePar(int aid)
{
  int startAtom=getRecordOffset();
  int endAtom=startAtom+getAtomsAssignedThisProc();
  if(!isAtomSent(aid))
  {
    if(aid>=startAtom && aid<endAtom)
      return molecule->atomvdwtype(aid-recordOffset);
  }
  if(newRecvAtom(aid))
  {
    FullAtom *recvAtom=this->recvAtomList.begin();
    long long myIdx=newRecvAtomIdx(aid);
    return recvAtom[myIdx].vdwType;
  }
  NAMD_bug("There is a  problem in random_velocities_parallel");
}

// value of sigId for a given atom id
short ParallelIOMgr::atomsigIDPar(int aid)
{
  int startAtom=getRecordOffset();
  int endAtom=startAtom+getAtomsAssignedThisProc();
  if(!isAtomSent(aid))
  {
    if(aid>=startAtom && aid<endAtom)
      return molecule->getAtomSigId(aid-recordOffset);
  }
  if(newRecvAtom(aid))
  {
    FullAtom *recvAtom=this->recvAtomList.begin();
    long long myIdx=newRecvAtomIdx(aid);
    return recvAtom[myIdx].sigId;
  }
  NAMD_bug("There is a  problem in random_velocities_parallel");
}

// value of exclSigId for a given atom id
short ParallelIOMgr::atomexclsigIDPar(int aid)
{
  int startAtom=getRecordOffset();
  int endAtom=startAtom+getAtomsAssignedThisProc();
  if(!isAtomSent(aid))
  {
    if(aid>=startAtom && aid<endAtom)
      return molecule->getAtomExclSigId(aid-recordOffset);
  }
  if(newRecvAtom(aid))
  {
    FullAtom *recvAtom=this->recvAtomList.begin();
    long long myIdx=newRecvAtomIdx(aid);
    return recvAtom[myIdx].exclId;
  }
  NAMD_bug("There is a  problem in random_velocities_parallel");
}

// value of charge for a given atom id
unsigned int ParallelIOMgr::atomchargePar(int aid)
{
  int startAtom=getRecordOffset();
  int endAtom=startAtom+getAtomsAssignedThisProc();
  if(!isAtomSent(aid))
  {
    if(aid>=startAtom && aid<endAtom)
      return molecule->getEachAtomCharge(aid-recordOffset);
  }
  if(newRecvAtom(aid))
  {
    FullAtom *recvAtom=this->recvAtomList.begin();
    long long myIdx=newRecvAtomIdx(aid);
    return recvAtom[myIdx].charge;
  }
  NAMD_bug("There is a  problem in random_velocities_parallel");
}

// value of mass for a given atom id
unsigned int ParallelIOMgr::atommassPar(int aid)
{
  int startAtom=getRecordOffset();
  int endAtom=startAtom+getAtomsAssignedThisProc();
  bool found=false;
  if(!isAtomSent(aid))
  {
    if(aid>=startAtom && aid<endAtom)
    {
      found=true;
      return molecule->getEachAtomMass(aid-recordOffset);
    }

  }
  if(newRecvAtom(aid))
  {
    FullAtom *recvAtom=this->recvAtomList.begin();
    long long myIdx=newRecvAtomIdx(aid);
    return recvAtom[myIdx].mass;
  }
}

// helper method for calculating atom-to-patch assignment
void ParallelIOMgr::initAtomsPerPatchParallel()
{
  numProcsSentNumAtomsPerPatch=0;
  PatchMap *patchMap=Node::Object()->getPatchMap();
  atomsPerPatchParallel=new int[patchMap->numPatches()];
  for(int p=0;p<patchMap->numPatches();p++)
    atomsPerPatchParallel[p]=0;
  calcAtomsEachPatch();
}

// method for updating exclusions. Some options are not covered like alchFepOn etc
void ParallelIOMgr::updateExclusions()
{
  molecule->numTotalExclusions/=2;
  molecule->setNumCalcExclusions(molecule->getNumCalcExclusions()/2);
}

#endif

// accumulating the counters on pe0
void ParallelIOMgr::updateCounters(int numBonds,int numAngles, int numDihedrals, int numImpropers, int numCrossterms, int numTotalExclusions,int numCalcExclusions)
{
#ifdef MEM_OPT_VERSION
  molecule->numBonds+=numBonds;
  molecule->numDihedrals+=numDihedrals;
  molecule->numAngles+=numAngles;
  molecule->numImpropers+=numImpropers;
  molecule->numCrossterms+=numCrossterms;
  molecule->numTotalExclusions+=numTotalExclusions;
  molecule->setNumCalcExclusions(molecule->getNumCalcExclusions()+numCalcExclusions);
#endif
}


#ifdef MEM_OPT_VERSION
// updating hydrogen group info after migration group exchange
void ParallelIOMgr::updateHydrogenGP()
{
  newTotalAtoms=this->getAtomsAssignedThisProc()-atomsLeavingProc+atomsEnteringProc;
  hydrogenGroupPar.resize(newTotalAtoms);
  HydrogenGroupID *hgPar = hydrogenGroupPar.begin();
  HydrogenGroupID *hg = molecule->hydrogenGroup.begin();
  int count=0;
  for(int i=0;i<atomsAssignedThisProc;i++)
  {
    int atomID=hg[i].atomID;
    int startAtom=this->getRecordOffset();
    int endAtom=startAtom+this->getAtomsAssignedThisProc();
    // this atom was not sent and still part of this input proc
    if(!isAtomSent(atomID))
    {
      if(atomID>=startAtom && atomID<endAtom)
      {
	hgPar[count].atomID=hg[i].atomID;  // currently unsorted
	hgPar[count].atomsInGroup=hg[i].atomsInGroup;  // currently only 1 in group                
        hgPar[count].isGP=hg[i].isGP;  // assume it is a group parent
        hgPar[count].GPID=hg[i].GPID;  // assume it is a group parent
        hgPar[count].sortVal=hg[i].sortVal;  // for group sorting
	hgPar[count].isMP=hg[i].isMP;
	hgPar[count].atomsInMigrationGroup=hg[i].atomsInMigrationGroup;
	hgPar[count].MPID=hg[i].MPID;
	hgPar[count].waterVal=hg[i].waterVal;
        count++;
      }
    }
  }
  FullAtom *recvAtom=this->recvAtomList.begin();
  for(int i=0;i<atomsEnteringProc;i++)
  {
    hgPar[count].atomID=recvAtom[i].id;  // currently unsorted
    hgPar[count].atomsInGroup=recvAtom[i].atomsInGroup;  // currently only 1 in group
    hgPar[count].GPID=recvAtom[i].GPID;  // assume it is a group parent
    hgPar[count].sortVal=recvAtom[i].sortVal;  // for group sorting
    hgPar[count].isGP = 1;
    hgPar[count].isMP=recvAtom[i].isMP;
    hgPar[count].atomsInMigrationGroup=recvAtom[i].migrationGroupSize;
    hgPar[count].MPID=recvAtom[i].MPID;
    hgPar[count].waterVal=recvAtom[i].waterVal;
    if(hgPar[count].atomsInGroup==0) hgPar[count].isGP = 0;
    count++;
  }
}

// helper method for calling method which updates hydrogen group structure after migration exchange
void ParallelIOMgr::updateHydrogenGroup()
{
  updateHydrogenGP();
  molecule->numCalcExclusions=molecule->numTotalExclusions;
  CProxy_ParallelIOMgr pIO(thisgroup);
  if(CkMyPe())
    pIO[0].updateCounters(molecule->numBonds,molecule->numAngles,molecule->numDihedrals,molecule->numImpropers,molecule->numCrossterms,molecule->numTotalExclusions,molecule->numCalcExclusions);
}

// method for reading molecule file
void ParallelIOMgr::readMolecule()
{
  molecule->read_compressed_psf_file_parallelIO((state->getConfigList()->find("structure"))->data,Node::Object()->parameters);
}
#endif

ParallelIOMgr::ParallelIOMgr()
{
#ifdef MEM_OPT_VERSION
  CkpvAccess(BOCclass_group).ioMgr = thisgroup;
  numInputProcs=-99;
  procsReceived=0;
#endif
}
#ifdef MEM_OPT_VERSION
// init method which is called from scripttcl.C
void ParallelIOMgr::initParallelIOMgr(long prc,long numinputprocs,long numRecords,char *pdb,bool flip,int *inputProcArr){
  procsReceived=0;
  CkpvAccess(BOCclass_group).ioMgr = thisgroup;
  atomsLeavingProc=atomsEnteringProc=0;
  proc=prc;
  inputProcArray=inputProcArr;
  doFlip=flip;
  pdbFileName=pdb;
  totalAtomsSys=numRecords;
  numInputProcs=numinputprocs;

	// remainderAtoms will have the remaining atoms if the total num of atoms are not divisible with the num of output procs. These aditional atoms
	// are assigned to the last processor.
  long remainderAtoms=totalAtomsSys%numInputProcs;
  long atomsPerProc=totalAtomsSys/numInputProcs;
	// this atomsPerProcNew is different only is the total number of atoms are not divisible by the number of output procs.
  atomsInLastProc=atomsPerProc;
  atomsPerProcOnetoNMinus1=atomsPerProc;
  atomsAssignedThisProc=atomsPerProc;


	// this to assign the remaining atoms to the last proc
  if(remainderAtoms>0) atomsInLastProc+=remainderAtoms; 
  if(proc==inputProcArray[numinputprocs-1]) atomsAssignedThisProc=atomsInLastProc;
	
  int inputProcIdx=-1;
  for(int i=0;i<numInputProcs;i++)
  {
    if(inputProcArray[i]==proc)
    {
      inputProcIdx=i;
      break;
    }
  }
  recordOffset=atomsPerProcOnetoNMinus1*inputProcIdx;
}


// method called for doing the migration atom exchange
void ParallelIOMgr::sendAtomsToMigrationGpParents()
{
  hydrogenGPUpated=false;
  int *countArr;
  countArr=new int[numInputProcs];
  atomsLeavingProc=0;
  MoveFullAtomsMsg **msg=new MoveFullAtomsMsg*[numInputProcs];
  FullAtomList *migrationAtoms=new FullAtomList[numInputProcs];
  FullAtom *migList=migrationList->begin();
  FullAtom *a=migrationList->begin();
  int procsToSend=0;
  for(int j=0;j<numInputProcs;j++)
  {
    int count=0;
    int inputProcNum=inputProcArray[j];
    for(int i=0;i<getNumMigrationAtoms();i++)
    {
      if(migList[i].readProc==CkMyPe() && inputProcNum==migList[i].destProc)
      {
	FullAtom at;
        at.readProc=a[i].readProc;
        at.destProc=a[i].destProc;
        migrationAtoms[j].add(at);
        count++;
      }
    }
    countArr[j]=count;
    FullAtom *atomsThisProc=migrationAtoms[j].begin();

    if(countArr[j]>0)
    {
      HydrogenGroupID *hg = molecule->hydrogenGroup.begin();
      count=0;
      for(int i=0;i<getNumMigrationAtoms();i++)
      {
	if(migList[i].readProc==CkMyPe() && inputProcNum==migList[i].destProc)
        {
	  int atomIDIdx=migList[i].id%getAtomsPerProcOnetoNMinus1();
	  int atomDiv=migList[i].id%getAtomsPerProcOnetoNMinus1();
          if(atomDiv>=numInputProcs)
	    atomIDIdx+=getAtomsPerProcOnetoNMinus1();
        
	  CoordinateRecord *recs=getRecords();
	  atomsThisProc[count].id=migList[i].id;
          atomsThisProc[count].GPID=hg[migList[i].id-recordOffset].GPID;
	  atomsThisProc[count].isMP=hg[migList[i].id-recordOffset].isMP;
	  atomsThisProc[count].migrationGroupSize=hg[migList[i].id-recordOffset].atomsInMigrationGroup;
	  atomsThisProc[count].MPID=hg[migList[i].id-recordOffset].MPID;
	  atomsThisProc[count].waterVal=hg[migList[i].id-recordOffset].waterVal;
	  atomsThisProc[count].atomsInGroup=hg[migList[i].id-recordOffset].atomsInGroup;
          atomsThisProc[count].sortVal=hg[migList[i].id-recordOffset].sortVal;
          AtomCstInfo *ats=molecule->getAtoms();
          atomsThisProc[count].mass=molecule->getEachAtomMass(migList[i].id-recordOffset);
          atomsThisProc[count].charge=molecule->getEachAtomCharge(migList[i].id-recordOffset);
          atomsThisProc[count].sigId=molecule->getAtomSigId(migList[i].id-recordOffset);
          atomsThisProc[count].exclId=molecule->getAtomExclSigId(migList[i].id-recordOffset);
          atomsThisProc[count].vdwType=ats[migList[i].id-recordOffset].vdw_type;
          atomsThisProc[count].position.x=recs[atomIDIdx].x;
          atomsThisProc[count].position.y=recs[atomIDIdx].y;
          atomsThisProc[count].position.z=recs[atomIDIdx].z;

          if ( simParameters->initialTemp < 0.0 ) {
	    atomsThisProc[count].velocity.x=atomArrayParallelInputVel[atomIDIdx].x;
	    atomsThisProc[count].velocity.y=atomArrayParallelInputVel[atomIDIdx].y;
	    atomsThisProc[count].velocity.z=atomArrayParallelInputVel[atomIDIdx].z;
          }
          count++;
        }
      }
      msg[procsToSend] = new MoveFullAtomsMsg(inputProcNum,CkMyPe(), migrationAtoms[j]);
      atomsLeavingProc+=count;
      procsToSend++;
    }
  }
  atomsLeaving=new int[atomsLeavingProc];
  int count=0;
  for(int i=0;i<getNumMigrationAtoms();i++)
  {
    if(migList[i].readProc==CkMyPe())
    {
      atomsLeaving[count]=migList[i].id;
      count++;
    }
  }

  for(int i=0;i<procsToSend;i++)
  {
    CProxy_ParallelIOMgr pIO(thisgroup);
    pIO[msg[i]->toProc].receiveMigrationChildren((MoveFullAtomsMsg*)MoveFullAtomsMsg::pack(msg[i]));
  }

  if(checkForSendAndRecv())
    hydrogenGPUpated=true;
}

// helper method which calculates where each atom needs to move for migration atom exchange
void ParallelIOMgr::redistributionAtomInfoParallel2()
{
  AtomCstInfo *atoms=this->molecule->getAtoms();
  setNumMigrationAtoms(0);
  HydrogenGroupID *hg = molecule->hydrogenGroup.begin();
  for(int i=0;i<atomsAssignedThisProc;i++)
  {
    if(atoms[i].status==HydrogenAtom)
    {
      int parent=hg[i].MPID;
      int inputProcIdxForHead=parent/getAtomsPerProcOnetoNMinus1();
      int inputProcIdx=hg[i].atomID/getAtomsPerProcOnetoNMinus1();
      if(inputProcIdx==numInputProcs) inputProcIdx=inputProcIdx-1;
      if(inputProcIdxForHead==numInputProcs) inputProcIdxForHead=inputProcIdxForHead-1;
      if(inputProcIdxForHead!=inputProcIdx)
	setNumMigrationAtoms(getNumMigrationAtoms()+1);
    }
  }
  migrationList = new FullAtomList;
  for(int i=0;i<getNumMigrationAtoms();i++)
  {
    FullAtom a;
    a.readProc=-1;
    a.destProc=-1;
    migrationList->add(a);	      
  }

  int count=0;
  hg = molecule->hydrogenGroup.begin();
  FullAtom *migList=migrationList->begin();
  for(int i=0;i<atomsAssignedThisProc;i++)
  {
    if(atoms[i].status==HydrogenAtom)
    {
      int parent=hg[i].MPID;
      int inputProcIdxForHead=parent/getAtomsPerProcOnetoNMinus1();
      int inputProcIdx=hg[i].atomID/getAtomsPerProcOnetoNMinus1();
      if(inputProcIdx==numInputProcs) inputProcIdx=inputProcIdx-1;
      if(inputProcIdxForHead==numInputProcs) inputProcIdxForHead=inputProcIdxForHead-1;
      if(inputProcIdxForHead!=inputProcIdx)
      {
	migList[count].id=hg[i].atomID;
        migList[count].destProc=inputProcArray[inputProcIdxForHead];
        migList[count].readProc=inputProcArray[inputProcIdx];
        count++;
      }
    }
  }
}
            
// getting the position for atom recvd as a result of migration exchange
void ParallelIOMgr::get_position_for_atom_ParallelVel(Vector *pos, int aid){
  for(int i=0;i<atomsLeavingProc;i++)
  {
    if(aid==atomsLeaving[i])
      NAMD_bug("This atom is no longer part of this input proc (in PDB.C)");
  }
  FullAtom *recvAtom=this->recvAtomList.begin();
  long long myIdx=newRecvAtomIdx(aid);
  if(myIdx!=-1)
  {
    pos->x=recvAtom[myIdx].velocity.x;
    pos->y=recvAtom[myIdx].velocity.y;
    pos->z=recvAtom[myIdx].velocity.z;
    return;	    
  }

  int atomNo=aid%atomsPerProcOnetoNMinus1;
  int atomDiv=aid/atomsPerProcOnetoNMinus1;
  if(atomDiv>=numInputProcs)
    atomNo+=atomsPerProcOnetoNMinus1;

  pos->x = atomArrayParallelInputVel[atomNo].x;
  pos->y = atomArrayParallelInputVel[atomNo].y;
  pos->z = atomArrayParallelInputVel[atomNo].z;

}

// method used to read velocity file
void ParallelIOMgr::readAtomInfoParallelInputVel(const char *pdbfilename,int proc)
{
  atomArrayParallelInputVel = new Velocity[atomsAssignedThisProc];
  if ( atomArrayParallelInputVel == NULL )
    NAMD_die("memory allocation failed in PDB::PDB");
  FILE *fp;
  if ( (fp = fopen(pdbfilename, "rb")) == NULL)
  {
    char errmsg[256];
     printf(errmsg, "Uggggggggggggggggnable 1111111  to open binary file 111111111 name=/u/ac/sarood/apoa1/osman.dcd\n" );
  }

  off_t startbyte=sizeof(int)+recordOffset*24;        
  for(int i=0;i<atomsAssignedThisProc;i++)
  {
    fseek(fp,startbyte+i*24,SEEK_SET);
    BigReal filen[3];
    if (fread(filen, sizeof(double), 3, fp) != (size_t)3)
    {
      char errmsg[256];
      printf("Error reading file in ParallelIOMgr\n");
    }
    else
    {
      if(doFlip) flipNum((char *)filen, sizeof(double), 3);
      atomArrayParallelInputVel[i].x=filen[0];
      atomArrayParallelInputVel[i].y=filen[1];
      atomArrayParallelInputVel[i].z=filen[2];
    }

  }
}

// delete storage created for final atom exchange
void ParallelIOMgr::deleteParallelStorage(bool isVel)
{
  delete [] records;
  if(atomsLeavingProc>0)  delete [] atomsLeaving;
}

// get atom position from parallel IO mgr
void ParallelIOMgr::get_position_for_atom_Parallel(Vector *pos, int aid){
  for(int i=0;i<atomsLeavingProc;i++)
  {
    if(aid==atomsLeaving[i])
      NAMD_bug("This atom is no longer part of this input proc (in PDB.C)");
  }
  FullAtom* recvAtom=this->recvAtomList.begin();
  long long myIdx=newRecvAtomIdx(aid);
  if(myIdx!=-1) 
  {
    pos->x=recvAtom[myIdx].position.x;
    pos->y=recvAtom[myIdx].position.y;
    pos->z=recvAtom[myIdx].position.z;
    return;
  }
  int atomNo=aid%atomsPerProcOnetoNMinus1;
  int atomDiv=aid/atomsPerProcOnetoNMinus1;
  if(atomDiv>=numInputProcs)
    atomNo+=atomsPerProcOnetoNMinus1;
  pos->x = records[atomNo].x;
  pos->y = records[atomNo].y;
  pos->z = records[atomNo].z;
}

//method used to read position
void ParallelIOMgr::readCoordinates()
{
  records = new CoordinateRecord[atomsAssignedThisProc];
  if ( records == NULL )
  {
    NAMD_die("memory allocation failed in PDB::PDB");
  }

  FILE *fp;    //  File descriptor
  char *coorFileN;
  coorFileN=pdbFileName;

  fp = fopen(coorFileN, "rb");

  if (! fp) {
    char s[500];
    sprintf(s, "The coordinate restart file '%s' specified as parameter 'bincoordinates' in config file can not be opened", pdbFileName);
    NAMD_err(s);
  }
  long long startbyte1=sizeof(int)+recordOffset*24;
  off_t startbyte=sizeof(int)+recordOffset*24;
  double t0=CmiWallTimer();
  fseek(fp,startbyte,SEEK_SET);//77 is the length of one line in bytes
  double t1=CmiWallTimer();	
  fflush(stdout);
  for(int i=0;i<atomsAssignedThisProc;i++)
  {
//    fseek(fp,startbyte+i*24,SEEK_SET);//77 is the length of one line in bytes
    BigReal filen[3];
    if (fread(filen, sizeof(double), 3, fp) != (size_t)3)
    {
      char errmsg[256];
      printf(errmsg, "Error reading binary file 22222222 11111 aaaaaaaaaaaa");
    }
    else
    {
      if(doFlip) flipNum((char *)filen, sizeof(double), 3);
      records[i].x=filen[0];
      records[i].y=filen[1];
      records[i].z=filen[2];
    }
  }
  double t2=CmiWallTimer();
//printf("+++++ PROC#%d seek time=%f read time=%f total time=%f input file=%s ++++++++++\n",proc,t1-t0,t2-t1,t2-t0,coorFileN);
}
#endif
ParallelIOMgr::~ParallelIOMgr(void){

}
// Packs the linked list associated with ConfigList
void* ConfigListMessage::pack(ConfigListMessage* inmsg)
{
  char *buf=inmsg->cfgList->packHelper(inmsg);
  delete inmsg;
  return (void*) buf;
}


//To unpack the ConfigList at each input proc
ConfigListMessage* ConfigListMessage::unpack(void *inbuf)
{

  char *buf=(char*)inbuf;
  ConfigListMessage *pmsg=(ConfigListMessage*)CkAllocBuffer(inbuf, sizeof(ConfigListMessage));
  pmsg->cfgList=new ConfigList();
  ConfigList *cfgLst=pmsg->cfgList;

  int readStart=0;char test[100];
  int numParams,atoms;
  bool flip;
  memcpy(&numParams,buf,sizeof(int));
  memcpy(&atoms,(buf+sizeof(int)),sizeof(int));
  memcpy(&flip,(buf+sizeof(int)*2),sizeof(bool));
  pmsg->doFlip=flip;
  readStart+=sizeof(int)*2+sizeof(bool);
  pmsg->numAtoms=atoms;
  for(int i=0;i<numParams;i++)
  {
    int nameLen;char name[100];
    memcpy(&nameLen,(void*)(buf+readStart),sizeof(int));
    readStart+=sizeof(int);
    memcpy(name,(void*)(buf+readStart),nameLen);
    name[nameLen]='\0';
    readStart+=nameLen;

    int numData;
    memcpy(&numData,(void*)(buf+readStart),sizeof(int));
    readStart+=sizeof(int);
    for(int j=0;j<numData;j++)
    {
      int dataLen;char data[100];
      memcpy(&dataLen,(void*)(buf+readStart),sizeof(int));
      readStart+=sizeof(int);
      memcpy(data,(void*)(buf+readStart),dataLen);
      data[dataLen]='\0';
      readStart+=dataLen;
      cfgLst->add_element(name,nameLen,data,dataLen);
    }
  }
  CkFreeMsg(inbuf);
  return pmsg;
}
PACK_MSG(AtomsPerPatchMsg,
  PACK_AND_NEW_ARRAY(patchAtomCnt,numPatches);
  PACK(numPatches);
)
PACK_MSG(NodeAssignmentMsg,
  PACK_AND_NEW_ARRAY(nodeAssigned,totalPatches);
  PACK_AND_NEW_ARRAY(atomsPerPatchParallel,totalPatches);
  PACK(totalPatches);
)

PACK_MSG(MoveFullAtomsMsg,
  PACK(toProc);
  PACK(fromProc);
  PACK_RESIZE(atom);
)

PACK_MSG(AtomListForPatchCreateMsg,
  PACK_RESIZE(atom);
  PACK(patchesToSend);
  PACK_AND_NEW_ARRAY(atomsPerPatch,patchesToSend);
  PACK_AND_NEW_ARRAY(patcheNumsToSend,patchesToSend);
)


#include "ParallelIOMgr.def.h"
