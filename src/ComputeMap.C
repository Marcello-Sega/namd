/**
***  Copyright (c) 1995, 1996, 1997, 1998, 1999, 2000 by
***  The Board of Trustees of the University of Illinois.
***  All rights reserved.
**/

#include <stdlib.h>
#include <stdio.h>

#include "charm++.h"

#include "ComputeMap.h"
#include "Compute.h"

#define MIN_DEBUG_LEVEL 4
//#define DEBUGM
#include "Debug.h"

// Singleton method
ComputeMap *ComputeMap::Instance() {
  if (CpvAccess(ComputeMap_instance) == 0) {
    CpvAccess(ComputeMap_instance) = new ComputeMap;	// this is never deleted
  }
  return CpvAccess(ComputeMap_instance);
}


//----------------------------------------------------------------------
ComputeMap::ComputeMap(void)
{
  nComputes=0;
  nPatchBased=0;
  nAtomBased=0;
  nAllocated=0;
  computeData = NULL;
}

//----------------------------------------------------------------------
ComputeMap::~ComputeMap(void)
{
  if (computeData != NULL)
  {
    for(int i=0; i<nComputes; i++)
    {
	delete[] computeData[i].pids;	// alloced and realloced above
    }
    delete[] computeData;	// alloced and realloced above
  }    
}

void
ComputeMap::checkMap(void)
{
  int computeCount = nComputes;
  for (int i=0; i<nComputes; i++) {
    if (computeData[i].compute) {
      computeCount++;
      if (! (computeData[i].compute->cid == i)) {
	DebugM(4, "ComputeID("<<computeData[i].compute->cid<<") != ComputeID("
	  << i <<")\n");
      }
    }
  }
  DebugM(4, "Compute Count = " << computeCount << "\n");
}

#undef PACK
#define PACK(type,data) { memcpy(b, &data, sizeof(type)); b += sizeof(type); }

int ComputeMap::packSize(void)
{
  int i;
  int size = 0;
  size += 4 * sizeof(int);
  for(i=0;i<nComputes;++i)
  {
    size += sizeof(ComputeData);
    size += computeData[i].numPidsAllocated * sizeof(PatchRec);
  }
  return size;
}

void ComputeMap::pack (char *buffer)
{
  DebugM(4,"Packing ComputeMap\n");
  int i,j;

  // fill in the data
  char *b = buffer;
  PACK(int,nPatchBased);
  PACK(int,nAtomBased);
  PACK(int,nComputes);
  PACK(int,nAllocated);
  for(i=0;i<nComputes;++i)
  {
    PACK(ComputeData,computeData[i]);
    for(j=0;j<computeData[i].numPidsAllocated;++j)
      PACK(PatchRec,computeData[i].pids[j]);
  }
}

#undef UNPACK
#define UNPACK(type,data) { memcpy(&data, b, sizeof(type)); b += sizeof(type); }

void ComputeMap::unpack (char *ptr)
{
  // Must copy over the Compute * to new ComputeMap! 
  ComputeData *oldComputeData = computeData;
  int oldNComputes = nComputes;
  computeData = NULL;

  DebugM(4,"Unpacking ComputeMap\n");
  int i,j;
  char *b = (char*)ptr;
  UNPACK(int,nPatchBased);
  UNPACK(int,nAtomBased);
  UNPACK(int,nComputes);
  UNPACK(int,nAllocated);
  computeData = new ComputeData[nAllocated];	// deleted during ~.
  for(i=0;i<nComputes;++i)
  {
    UNPACK(ComputeData,computeData[i]);
    computeData[i].pids = new PatchRec[computeData[i].numPidsAllocated]; // deleted during ~
    for(j=0;j<computeData[i].numPidsAllocated;++j)
      UNPACK(PatchRec,computeData[i].pids[j]);
  }

  if (oldComputeData) {
    if (nComputes != oldNComputes) {
      NAMD_die("number of computes in new patchmap has changed!\n");
      return;
    }

    for (int i=0; i<nComputes; i++) {
      computeData[i].compute = oldComputeData[i].compute;
      delete[] oldComputeData[i].pids;
    }
    delete[] oldComputeData;
    oldComputeData=NULL;
  }
  DebugM(4,"Done Unpacking ComputeMap\n");
}


//----------------------------------------------------------------------
int ComputeMap::numComputes(void)
{
  return nComputes;
}

//----------------------------------------------------------------------
int ComputeMap::numPatchBased(void)
{
  return nPatchBased;
}

//----------------------------------------------------------------------
int ComputeMap::numAtomBased(void)
{
  return nAtomBased;
}

//----------------------------------------------------------------------
int ComputeMap::isPatchBased(ComputeID cid)
{
  if (computeData != NULL)
    return computeData[cid].patchBased;
  else return -1;
}

//----------------------------------------------------------------------
int ComputeMap::isAtomBased(ComputeID cid)
{
  if (computeData != NULL)
    return !computeData[cid].patchBased;
  else return -1;
}

//----------------------------------------------------------------------
int ComputeMap::node(ComputeID cid)
{
  if (computeData != NULL)
    return computeData[cid].node;
  else return -1;
}

void ComputeMap::setNode(ComputeID cid, NodeID node) {
  computeData[cid].node = node;
}

NodeID ComputeMap::newNode(ComputeID cid)
{
  return (computeData[cid].moveToNode);
}


void ComputeMap::setNewNode(ComputeID cid, NodeID node) {
  computeData[cid].moveToNode = node;
}

//----------------------------------------------------------------------
int ComputeMap::numPids(ComputeID cid)
{
  if (computeData != NULL)
    return computeData[cid].numPids;
  else return -1;
}

//----------------------------------------------------------------------
int ComputeMap::pid(ComputeID cid,int i)
{
  if ((computeData != NULL) && (i < computeData[cid].numPids))
    return computeData[cid].pids[i].pid;
  else return -1;
}

//----------------------------------------------------------------------
ComputeType ComputeMap::type(ComputeID cid)
{
  if (computeData != NULL)
    return computeData[cid].type;
  else return computeErrorType;
}

//----------------------------------------------------------------------
int ComputeMap::partition(ComputeID cid)
{
  if (computeData != NULL)
    return computeData[cid].partition;
  else return computeErrorType;
}
//----------------------------------------------------------------------
int ComputeMap::numPartitions(ComputeID cid)
{
  if (computeData != NULL)
    return computeData[cid].numPartitions;
  else return computeErrorType;
}

//----------------------------------------------------------------------
int ComputeMap::allocateCids(int n)
{
  if (computeData != NULL)
  {
    register int i;
    for(i=0; i<nComputes; i++)
    {
      if (computeData[i].pids != NULL)
      {
	delete [] computeData[i].pids;	// alloc in storeCompute and unpack
	computeData[i].pids=0;
      }
      delete [] computeData;	// delete before realloc

      computeData = NULL;
    }
  }
  nComputes = nPatchBased = nAtomBased = 0;

  // Constructor zero's out array elements
  computeData = new ComputeData[n];	// realloced after delete
  nAllocated = n;

  return 0;
}

//----------------------------------------------------------------------
ComputeID ComputeMap::storeCompute(int inode, int maxPids, 
				   ComputeType type, 
				   int partition,int numPartitions)
{
  int cid;

  if (!computeData)
    NAMD_die("ID's not allocated in ComputeMap::storeCompute().\n");
    //return -1;                   // Have to allocate first

  if (nComputes == nAllocated)
    NAMD_die("Allocated ID's exceeded in ComputeMap::storeCompute().\n");
    //return -1;                   // Used up all the allocated entries

  cid = nComputes;
  nComputes++;

  computeData[cid].node=inode;

  computeData[cid].type = type;
  computeData[cid].partition = partition;
  computeData[cid].numPartitions = numPartitions;

  computeData[cid].patchBased = true;
  nPatchBased++;

  computeData[cid].numPids = 0;
  computeData[cid].pids = new PatchRec[maxPids]; // delete in allocCids and ~
  computeData[cid].numPidsAllocated = maxPids;

  return cid;
}

//----------------------------------------------------------------------
int ComputeMap::newPid(ComputeID cid, PatchID pid, int trans)
{
  if (!computeData)
    return -1;                   // Have to allocate first

  if ((cid < 0) || (cid >= nComputes))
    return -1;                   // Have to store the cid first

  if (computeData[cid].numPids == computeData[cid].numPidsAllocated)
    return -1;                   // Out of space for dependents

  computeData[cid].pids[computeData[cid].numPids].pid=pid;
  computeData[cid].pids[computeData[cid].numPids].trans=trans;
  computeData[cid].numPids++;
  return 0;
}

//----------------------------------------------------------------------
void ComputeMap::printComputeMap(void)
{
  DebugM(2,"---------------------------------------");
  DebugM(2,"---------------------------------------\n");

  DebugM(2,"nComputes = " << nComputes << '\n');
  DebugM(2,"nPatchBased = " << nPatchBased << '\n');
  DebugM(2,"nAtomBased = " << nAtomBased << '\n');
  DebugM(2,"nAllocated = " << nComputes << '\n');
  for(int i=0; i < nComputes; i++)
  {
    DebugM(2,"Compute " << i << '\n');
    DebugM(2,"  node = " << computeData[i].node << '\n');
    DebugM(2,"  patchBased = " << computeData[i].patchBased << '\n');
    DebugM(2,"  numPids = " << computeData[i].numPids << '\n');
    DebugM(2,"  numPidsAllocated = " << computeData[i].numPidsAllocated << '\n');
    for(int j=0; j < computeData[i].numPids; j++)
    {
      DebugM(2,computeData[i].pids[j].pid);
      if (!((j+1) % 6))
	DebugM(2,'\n');
    }
    DebugM(2,"\n---------------------------------------");
    DebugM(2,"---------------------------------------\n");

  }
}

