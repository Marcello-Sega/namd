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
class ComputeMapDistribMsg;

class WorkDistrib : public BOCclass
{
public:
  WorkDistrib();
  ~WorkDistrib(void);

  // static void messageMovePatchDone();
  // void movePatchDone();

  static void messageEnqueueWork(Compute *);
  void enqueueWork(LocalWorkMsg *msg);
  void enqueueBonds(LocalWorkMsg *msg);
  void enqueueAngles(LocalWorkMsg *msg);
  void enqueueDihedrals(LocalWorkMsg *msg);
  void enqueueImpropers(LocalWorkMsg *msg);
  void enqueuePme(LocalWorkMsg *msg);
  void enqueueSelfA(LocalWorkMsg *msg);
  void enqueueSelfB(LocalWorkMsg *msg);
  void enqueueWorkA(LocalWorkMsg *msg);
  void enqueueWorkB(LocalWorkMsg *msg);
  void enqueueWorkC(LocalWorkMsg *msg);

  void mapComputes(void);
  void sendMaps(void);
  void saveComputeMapChanges(int,int);
  void recvComputeMapChanges(ComputeMapChangeMsg *);
  void doneSaveComputeMap();
  void createHomePatches(void);
  void distributeHomePatches(void);
  void patchMapInit(void);
  void assignNodeToPatch(void);

  void saveMaps(MapDistribMsg *msg);

private:
  void mapComputeNonbonded(void);
  void mapComputeHomePatches(ComputeType);
  void mapComputePatch(ComputeType);
  void assignPatchesToLowestLoadNode(void);
  void assignPatchesRecursiveBisection(void);
  void assignPatchesRoundRobin(void);
  void velocities_from_PDB(char *filename, 
			   Vector *v, int totalAtoms);
  void velocities_from_binfile(char *fname, Vector *vels, int n);
  void random_velocities(BigReal Temp, Molecule *structure,
			 Vector *v, int totalAtoms);
  void remove_com_motion(Vector *vel, Molecule *structure, int n);

  Boolean mapsArrived;
  Boolean awaitingMaps;
  CthThread awaitingMapsTh;

  int saveComputeMapReturnEP;
  int saveComputeMapReturnChareID;
  int saveComputeMapCount;
};

#include <string.h>
#include "PatchMap.h"
#include "ComputeMap.h"

class ComputeMapDistribMsg : public CMessage_ComputeMapDistribMsg
{
public:
  ComputeMapDistribMsg(void) : computeMap(0) { ; }
  ComputeMap *computeMap;

  // pack and unpack functions
  static void* pack(ComputeMapDistribMsg* msg)
  {
    int computeMapSize = msg->computeMap->packSize();
    char *buffer = (char*)CkAllocBuffer(msg,computeMapSize);
    msg->computeMap->pack(buffer);
    delete msg;
    return buffer;
  }

  static ComputeMapDistribMsg* unpack(void *ptr)
  {
    void *_ptr = CkAllocBuffer(ptr, sizeof(ComputeMapDistribMsg));
    ComputeMapDistribMsg *m = new (_ptr) ComputeMapDistribMsg();
    m->computeMap = ComputeMap::Object();
    m->computeMap->unpack((char*)ptr);
    CkFreeMsg(ptr);
    return m;
  }
};

class MapDistribMsg : public CMessage_MapDistribMsg
{
public:
  MapDistribMsg(void) : patchMap(0), computeMap(0) { ; }
  PatchMap *patchMap;
  ComputeMap *computeMap;

  // pack and unpack functions
  static void* pack(MapDistribMsg *msg)
  {
    int patchMapSize = msg->patchMap->packSize();
    int computeMapSize = msg->computeMap->packSize();
    int length = sizeof(int) + patchMapSize + sizeof(int) + computeMapSize;
    char *buffer = (char*)CkAllocBuffer(msg,length);
    char *b = buffer;
    *((int*)b) = patchMapSize;
    b += sizeof(int);
    msg->patchMap->pack(b);
    b += patchMapSize;
    *((int*)b) = computeMapSize;
    b += sizeof(int);
    msg->computeMap->pack(b);
    b += computeMapSize;
    delete msg;
    return buffer;
  }

  static MapDistribMsg* unpack(void *ptr)
  {
    void *_ptr = CkAllocBuffer(ptr, sizeof(MapDistribMsg));
    MapDistribMsg *m = new (_ptr) MapDistribMsg;
    char *buffer = (char*)ptr;
    int patchMapSize = *((int*)buffer);
    buffer += sizeof(int);
    m->patchMap = PatchMap::Object();
    if ( ! ( m->patchMap->patchData ) )
    {
      m->patchMap->unpack(buffer);
    }
    buffer += patchMapSize;
    int computeMapSize = *((int*)buffer);
    buffer += sizeof(int);
    m->computeMap = ComputeMap::Object();
    if ( ! ( m->computeMap->computeData ) )
    {
      m->computeMap->unpack(buffer);
    }
    buffer += computeMapSize;
    CkFreeMsg(ptr);
    return m;
  }
};

class ComputeMapChangeMsg : public CMessage_ComputeMapChangeMsg
{
public:
  int newNodes[20000];
  int numNewNodes;
};


#endif /* WORKDISTRIB_H */

