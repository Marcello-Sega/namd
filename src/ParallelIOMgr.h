#ifndef PARALLELIOMGR_H
#define PARALLELIOMGR_H  


#include "ProcessorPrivate.h"
#include "charm++.h"
#include "BOCgroup.h"
#include "common.h"
#include "CompressPsf.h"
#include "Hydrogen.h"
#include "Vector.h"
#include "NamdState.h"
#include "ParallelIOMgr.decl.h"
#include "Node.decl.h"
#include "PatchMgr.h"

class CoordinateRecord;
class ConfigList;
//class Molecule;
class PatchMap;

class AtomListForPatchCreateMsg :public CMessage_AtomListForPatchCreateMsg{
public:
  long long patchesToSend;
  FullAtomList atom;
  int *atomsPerPatch;
  int *patcheNumsToSend;

    static void* pack(AtomListForPatchCreateMsg *msg);
    static AtomListForPatchCreateMsg* unpack(void *ptr);
};

class MoveFullAtomsMsg : public CMessage_MoveFullAtomsMsg {
public:
    int toProc;
    int fromProc;
    FullAtomList atom;

    MoveFullAtomsMsg(void) { ; }

    MoveFullAtomsMsg(int t,int f, FullAtomList a) :atom(a)
    {
      toProc=t;
      fromProc=f;
    }

  // pack and unpack functions
    static void* pack(MoveFullAtomsMsg *msg);
    static MoveFullAtomsMsg* unpack(void *ptr);
};

class NodeAssignmentMsg: public CMessage_NodeAssignmentMsg {
public:
  int *nodeAssigned;
  int *atomsPerPatchParallel;
  int totalPatches;
  static void* pack(NodeAssignmentMsg* msg);
  static NodeAssignmentMsg* unpack(void *ptr);
};
class AtomsPerPatchMsg: public CMessage_AtomsPerPatchMsg {
public:
  int *patchAtomCnt;
  int numPatches;
  static void* pack(AtomsPerPatchMsg* msg);
  static AtomsPerPatchMsg* unpack(void *ptr);
};
class ConfigListMessage : public CMessage_ConfigListMessage{

public:
  ConfigList *cfgList;
  int numAtoms;
  bool doFlip;
  static void* pack(ConfigListMessage *);
  static ConfigListMessage* unpack(void*);
};

class ParallelIOMgr : public BOCclass
{
private:
  int proc;
  long atomsPerProcOnetoNMinus1;
  long atomsInLastProc;
  long atomsAssignedThisProc;
  long numInputProcs;
  long numRecords;
  char *pdbFileName;
  long long recordOffset;
  CoordinateRecord *records;
  bool doFlip;
  int *inputProcArray;
  int numMigrationAtoms;
  Molecule *molecule;
  NamdState *state;
  int *atomsPerPatchParallel;
  PatchMgr *patchMgr;
  int newTotalAtoms;
  SimParameters *simParameters;
  int procsReceived;
  int currentAtomsAfterPatchExchange;
  int *assignedNodesParallel;
  int numProcsSentNumAtomsPerPatch;
  int *patchesAfterPatchExchange; //these are populated when we get coorindates for patches assigned to this proc by other procs.    

public:

  void checkConfigParas();
  void atomExchangeForPatchCreation();
  void caclNumAtomsInEachPatchParallel();
  void receiveMigrationChildren(MoveFullAtomsMsg *msg);
  void receiveNumAtomsPerPatch(AtomsPerPatchMsg *msg1);
  void sendCoordinatesToProcs();
  void readCoordinatesAndVel(ConfigListMessage *pmsg);
  ParallelIOMgr();
  void getPatchCoordinates(AtomListForPatchCreateMsg *msg);
  void updateCounters(int numBonds,int numAngles, int numDihedrals, int numImpropers, int numCrossterms, int numTotalExclusions,int numCalcExclusions);
  void sendAssignedNode(NodeAssignmentMsg *msg);
#ifdef MEM_OPT_VERSION
  void updateExclusions();
  FullAtomList *migrationList;
  FullAtomList recvAtomList;
  FullAtomList finalAtomList;
  int *mappingHydGPPerProc;
  HydrogenGroup hydrogenGroupPerProc;
  bool hydrogenGPUpated;
  HydrogenGroup hydrogenGroupPar;
  int *hydrogenGroupParIdx;
  
  vector<int> *patchAtomList;
  vector<int> *patchAtomIdMappingList;
  void initAssignedNodesParallel();
  SimParameters *getSimParameters() {return simParameters;}
  void setDoFlip(bool s){doFlip=s;}
  int getNumInputProcs(){return numInputProcs;}
  int getNewTotalAtoms(){return newTotalAtoms;}
  NamdState *getState(){return state;}
  CProxy_Node cm;
  Bool is_hydrogenGroupParentPar(int anum)
  {
    return (hydrogenGroupPerProc[mappingHydGPPerProc[anum]].isGP);
  }

  int getIsMP(int anum)
  {
    return (hydrogenGroupPerProc[mappingHydGPPerProc[anum]].isMP);
  }

  int get_groupSizePar(int anum)
  {
    return (hydrogenGroupPerProc[mappingHydGPPerProc[anum]].atomsInGroup);
  }
  int get_MigrationGpSizePar(int anum)
  {
    return (hydrogenGroupPerProc[mappingHydGPPerProc[anum]].atomsInMigrationGroup);
  }
  int *getInputProcArray(){return inputProcArray;}
  void setInputProcArray(int *d){inputProcArray=d;}
  int totalAtomsSys;
  int inputProcIndex(int);
  void calcAtomsEachPatch();

  void initAtomsPerPatchParallel();
  short atomsigIDPar(int);
  short atomexclsigIDPar(int);
  short atomvdwtypePar(int aid);
  unsigned int atomchargePar(int aid);
  unsigned int atommassPar(int aid);
  int atomMigGrpSize(int);
  int atomIsMP(int);
  void mapHydGPPerProc();
  bool checkForSendAndRecv();
  bool isAtomSent(int);
  long long newRecvAtomIdx(int atomID);
  bool newRecvAtom(int);
  void initParallelInput();
  bool isInputProc(int);
  void initInputProcArray(int);
  void get_position_for_patch_creationVel(Vector *pos,int aid,int aidIdx);
  void get_position_for_patch_creation(Vector *pos,int aid,int aidIdx);
  void broadcastPatchInfo();
  void createAtomHoldPatchTransfer(int numPatches);
  void createPatchesThisProc();
  void deleteAtomHoldPatchTransfer();
  void updateHydrogenGP();
  void updateHydrogenGroup();
  void readMolecule();
  void initParallelIOMgr(long proc,long numInputProcs,long numRecords,char *pdb,bool flip,int *inputProcArr);
  int *getAtomsPerPatchParallel(){return atomsPerPatchParallel;}
  int getAtomsPerPatchParallel(int i){return atomsPerPatchParallel[i];}
  void setAtomsPerPatchParallel(int *p){atomsPerPatchParallel=p;}
  void incAtomsPerPatchParallel(int idx,int val){atomsPerPatchParallel[idx]=atomsPerPatchParallel[idx]+val;}
  Velocity *atomArrayParallelInputVel;
  void sendAtomsToMigrationGpParents();
  void get_position_for_atom_ParallelVel(Vector *pos, int aid);
  void deleteParallelStorage(bool isVel);
  int *atomsLeaving;
  int atomsLeavingProc,atomsEnteringProc;
  void readAtomInfoParallelInputVel(const char *pdbfilename,int proc);
  Molecule *getMolecule(){return molecule;}
  void setPointers(Molecule *m,CProxy_Node cm1,NamdState *st,PatchMgr *pm,SimParameters *s)
  {
    this->cm=cm1;this->molecule=m;this->state=st;this->patchMgr=pm;simParameters=s;
  }
  void redistributionAtomInfoParallel2();
  void setMolecule(Molecule *m){molecule=m;}
  int getNumMigrationAtoms(){return numMigrationAtoms;}
  void setNumMigrationAtoms(int i){numMigrationAtoms=i;}
  long long getRecordOffset() {return recordOffset;}
  void setRecordOffset(long long i) {recordOffset=i;}
  void get_position_for_atom_Parallel(Vector *pos, int aid);
  long getAtomsPerProcOnetoNMinus1()
  {
    return atomsPerProcOnetoNMinus1;
  }
  long getAtomsInLastProc()
  {
    return atomsInLastProc;
  }

  long getAtomsAssignedThisProc()
  {
    return atomsAssignedThisProc;
  }

  void setAtomsAssignedThisProc(int i)
  {
    atomsAssignedThisProc=i;
  }

  CoordinateRecord *getRecords()
  {
    return records;
  }

  void readCoordinates();
#endif
  ~ParallelIOMgr(void);
};

class CoordinateRecord{
public:
	BigReal x,y,z;
	CoordinateRecord()
	{
		x=y=z=0;
	}
	CoordinateRecord(BigReal x1,BigReal y1,BigReal z1)
	{
		x=x1;y=y1;z=z1;
	}

	void setCoor(BigReal x1,BigReal y1, BigReal z1)
	{
	  x=x1;y=y1;z=z1;
	}
};
#endif
