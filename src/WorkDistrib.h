/**
***  Copyright (c) 1995, 1996, 1997, 1998, 1999, 2000 by
***  The Board of Trustees of the University of Illinois.
***  All rights reserved.
**/

#ifndef _WORKDISTRIB_H
#define _WORKDISTRIB_H

#include "charm++.h"

#include "main.h"

#include "NamdTypes.h"
#include "BOCgroup.h"
#include "ComputeMap.h"
#include "WorkDistrib.decl.h"

class Node;
class Compute;
class Molecule;

// For Compute objects to enqueue themselves when ready to compute
class LocalWorkMsg : public CMessage_LocalWorkMsg
{
public:
  Compute *compute;
};

enum { maxPatchDepends = 126 };

class MapDistribMsg;
class ComputeMapChangeMsg;

class WorkDistrib : public BOCclass
{
public:
  WorkDistrib();
  ~WorkDistrib(void);

  // static void messageMovePatchDone();
  // void movePatchDone();

  static void messageEnqueueWork(Compute *);
  void enqueueWork(LocalWorkMsg *msg);
  void enqueueExcls(LocalWorkMsg *msg);
  void enqueueBonds(LocalWorkMsg *msg);
  void enqueueAngles(LocalWorkMsg *msg);
  void enqueueDihedrals(LocalWorkMsg *msg);
  void enqueueImpropers(LocalWorkMsg *msg);
  void enqueueThole(LocalWorkMsg *msg);  // Drude model
  void enqueueAniso(LocalWorkMsg *msg);  // Drude model
  void enqueueCrossterms(LocalWorkMsg *msg);
  void enqueuePme(LocalWorkMsg *msg);
  void enqueueSelfA(LocalWorkMsg *msg);
  void enqueueSelfB(LocalWorkMsg *msg);
  void enqueueWorkA(LocalWorkMsg *msg);
  void enqueueWorkB(LocalWorkMsg *msg);
  void enqueueWorkC(LocalWorkMsg *msg);
  void enqueueCUDA(LocalWorkMsg *msg);

  void mapComputes(void);
  void sendMaps(void);
  void saveComputeMapChanges(int,CkGroupID);
  void recvComputeMapChanges(ComputeMapChangeMsg *);
  void doneSaveComputeMap();

  FullAtomList *createAtomLists(void);
  void createHomePatches(void);

  int *caclNumAtomsInEachPatch(void);
  //create home patches without populating them with atoms' data
  //the only data set is the number atoms the patch contains
  void preCreateHomePatches(void);
  void fillOnePatchAtoms(int patchId, FullAtomList *onePatchAtoms, Vector *velocities);

  void distributeHomePatches(void);
  //populate the home patch with actual atom data
  void initAndSendHomePatch(void);

  void reinitAtoms(void);
  void patchMapInit(void);
  void assignNodeToPatch(void);

  void saveMaps(MapDistribMsg *msg);

private:
  void mapComputeNonbonded(void);
  void mapComputeNode(ComputeType);
  void mapComputeHomePatches(ComputeType);
  void mapComputeHomeTuples(ComputeType);
  void mapComputePatch(ComputeType);
  void assignPatchesToLowestLoadNode(void);
  void assignPatchesRecursiveBisection(void);
  void assignPatchesRoundRobin(void);
  void assignPatchesSpaceFillingCurve(void);
  void assignPatchesBitReversal(void);
  int  assignPatchesTopoGridRecBisection();

  void sortNodesAndAssign(int *assignedNode, int baseNodes = 0);
  void velocities_from_PDB(char *filename, 
			   Vector *v, int totalAtoms);
  void velocities_from_binfile(char *fname, Vector *vels, int n);
  void random_velocities(BigReal Temp, Molecule *structure,
			 Vector *v, int totalAtoms);
  void remove_com_motion(Vector *vel, Molecule *structure, int n);

  bool mapsArrived;
  bool awaitingMaps;
  CthThread awaitingMapsTh;

  int saveComputeMapReturnEP;
  CkGroupID saveComputeMapReturnChareID;
  int saveComputeMapCount;
};

#endif /* WORKDISTRIB_H */

