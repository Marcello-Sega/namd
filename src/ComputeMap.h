/**
***  Copyright (c) 1995, 1996, 1997, 1998, 1999, 2000 by
***  The Board of Trustees of the University of Illinois.
***  All rights reserved.
**/

#ifndef COMPUTEMAP_H
#define COMPUTEMAP_H

#include "NamdTypes.h"
#include "ProcessorPrivate.h"

class Compute;
class ComputeMgr;

enum ComputeType
{
  computeNonbondedSelfType,
  computeNonbondedPairType,
  computeBondsType,
  computeAnglesType,
  computeDihedralsType,
  computeImpropersType,
  computeSelfBondsType,
  computeSelfAnglesType,
  computeSelfDihedralsType,
  computeSelfImpropersType,
#ifdef DPMTA
  computeDPMTAType,
#endif
#ifdef DPME
  computeDPMEType,
#endif
  computePmeType,
  computeFullDirectType,
  computeGlobalType,
  computeEFieldType,
  computeSphericalBCType,
  computeCylindricalBCType,
  computeRestraintsType,
  computeErrorType
};

class ComputeMap
{
public:
  static ComputeMap *Instance();
  inline static ComputeMap *Object() { return CpvAccess(ComputeMap_instance); }

  void checkMap();

  ~ComputeMap(void);

  void registerCompute(ComputeID cid, Compute *c) {
    computeData[cid].compute = c;
    computeData[cid].moveToNode = -1;
  }

  // numComputes() returns the number of compute objects known
  // by the map.
  int numComputes(void);

  // numPatchBased() returns the number of compute objects
  // that are patch-based
  int numPatchBased(void);

  // numAtomBased() returns the number of compute objects
  // that are atom-based
  int numAtomBased(void);

  // isPatchBased(cid) returns true if the compute object
  // is patch based.
  int isPatchBased(ComputeID cid);

  // isAtomBased(cid) returns true if the compute object
  // is atom based.
  int isAtomBased(ComputeID cid);

  // node(cid) returns the node where the compute object currently exists.
  int node(ComputeID cid);

  void setNode(ComputeID cid, NodeID node);

  // newNode(cid,node) sets up map to tell WorkDistrib to send 
  // compute to new node
  NodeID newNode(ComputeID cid);

  void setNewNode(ComputeID cid, NodeID node);

  // numPids(cid) returns the number of patch ids which are registered
  // with this compute object.
  int numPids(ComputeID cid);
  
  // pid(cid,i) returns the i-th patch id registered
  // with the patch.  
  int pid(ComputeID cid, int i);

  // type(cid) returns the compute type of the given ComputeID
  ComputeType type(ComputeID cid);
  int partition(ComputeID cid);
  int numPartitions(ComputeID cid);

  // allocate_cids(n) tells the ComputeMap to set aside
  // room for n compute objects.  I need not use them all.
  int allocateCids(int n);

  // storeCompute(cid,node,maxPids) tells the ComputeMap to store
  // information about the indicated patch, and allocate space
  // for up to maxPids dependents
  ComputeID storeCompute(int node,int maxPids,ComputeType type,
			 int partition=0, int numPartitions=1);

  // newPid(cid,pid) stores the n patch ids associated with
  // compute id cid.
  int newPid(ComputeID cid, int pid, int trans = 13);

  void printComputeMap(void);

  Compute *compute(ComputeID cid) { return (computeData[cid].compute); };

  friend class ComputeMgr;

  struct PatchRec
  {
    PatchID pid;
    int trans;

    PatchRec() : pid(-1), trans(-1) { ; }
  };

  struct ComputeData
  {
    ComputeData() { 
      node = -1; moveToNode = -1; 
      patchBased = false; numPids = 0; numPidsAllocated = 0; 
      pids = NULL; compute = NULL; 
    }
    Compute *compute;
    int node;
    int moveToNode;
    ComputeType type;
    int partition;
    int numPartitions;
    Boolean patchBased;
    int numPids;
    int numPidsAllocated;
    PatchRec *pids;
  };
protected:
  friend class WorkDistrib;
  int packSize(void);
  void pack(char* buf);
  void unpack(char *buf);

  ComputeMap(void);

private:
  int nPatchBased;
  int nAtomBased;
  int nComputes;
  int nAllocated;
  ComputeData *computeData;
};

#endif /* COMPUTEMAP_H */

