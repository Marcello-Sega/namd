/***************************************************************************/
/*                                                                         */
/*              (C) Copyright 1996 The Board of Trustees of the            */
/*                          University of Illinois                         */
/*                           All Rights Reserved                           */
/*                                                                         */
/***************************************************************************/

/***************************************************************************
 * DESCRIPTION:
 *
 ***************************************************************************/

static char ident[] = "@(#)$Header: /home/cvs/namd/cvsroot/namd2/src/ComputeMap.C,v 1.1014 1997/11/07 20:17:36 milind Exp $";

#include <stdlib.h>
#include <stdio.h>

#include "chare.h"
#include "ckdefs.h"
#include "c++interface.h"

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

void * ComputeMap::pack (int *length)
{
  DebugM(4,"Packing ComputeMap\n");
  int i,j;

  // calculate memory needed
  int size = 0;
  size += 4 * sizeof(int);
  for(i=0;i<nComputes;++i)
  {
    size += sizeof(ComputeData);
    size += computeData[i].numPidsAllocated * sizeof(PatchRec);
  }
  *length = size;

  // allocate needed memory
  char * const buffer = new char[size];	// deleted in WorkDistrib.h::pack

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

  return buffer;
}

#undef UNPACK
#define UNPACK(type,data) { memcpy(&data, b, sizeof(type)); b += sizeof(type); }

void ComputeMap::unpack (void *in)
{
  // Must copy over the Compute * to new ComputeMap! 
  ComputeData *oldComputeData = computeData;
  int oldNComputes = nComputes;
  computeData = NULL;

  DebugM(4,"Unpacking ComputeMap\n");
  int i,j;
  char *b = (char*)in;
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
      iout << iPE << iERRORF 
        << "number of computes in new patchmap has changed!\n" << endi;
      CharmExit();
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
ComputeID ComputeMap::storeCompute(int inode, int maxPids, ComputeType type)
{
  int cid;

  if (!computeData)
    return -1;                   // Have to allocate first

  if (nComputes == nAllocated)
    return -1;                   // Used up all the allocated entries

  cid = nComputes;
  nComputes++;

  computeData[cid].node=inode;

  computeData[cid].type = type;

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
/***************************************************************************
 * RCS INFORMATION:
 *
 *	$RCSfile: ComputeMap.C,v $
 *	$Author: milind $	$Locker:  $		$State: Exp $
 *	$Revision: 1.1014 $	$Date: 1997/11/07 20:17:36 $
 *
 ***************************************************************************
 * REVISION HISTORY:
 *
 * $Log: ComputeMap.C,v $
 * Revision 1.1014  1997/11/07 20:17:36  milind
 * Made NAMD to run on shared memory machines.
 *
 * Revision 1.1013  1997/04/10 09:13:49  ari
 * Final debugging for compute migration / proxy creation for load balancing.
 * Lots of debug code added, mostly turned off now.
 * Fixed bug in PositionBox when Patch had no dependencies.
 * Eliminated use of cout and misuse of iout in numerous places.
 *                                            Ari & Jim
 *
 * Revision 1.1012  1997/04/08 07:08:12  ari
 * Modification for dynamic loadbalancing - moving computes
 * Still bug in new computes or usage of proxies/homepatches.
 * Works if ldbStrategy is none as before.
 *
 * Revision 1.1011  1997/04/04 23:34:17  milind
 * Got NAMD2 to run on Origin2000.
 * Included definitions of class static variables in C files.
 * Fixed alignment bugs by using memcpy instead of assignment in
 * pack and unpack.
 *
 * Revision 1.1010  1997/03/27 20:25:40  brunner
 * Changes for LdbCoordinator, the load balance control BOC
 *
 * Revision 1.1009  1997/03/20 23:53:38  ari
 * Some changes for comments. Copyright date additions.
 * Hooks for base level update of Compute objects from ComputeMap
 * by ComputeMgr.  Useful for new compute migration functionality.
 *
 * Revision 1.1008  1997/03/14 21:40:08  ari
 * Reorganized startup to make possible inital load
 * balancing by changing methods in WorkDistrib.
 * Also made startup more transparent and easier
 * to modify.
 *
 * Revision 1.1007  1997/03/04 22:37:07  ari
 * Clean up of code.  Debug statements removal, dead code removal.
 * Minor fixes, output fixes.
 * Commented some code from the top->down.  Mainly reworked Namd, Node, main.
 *
 * Revision 1.1006  1997/02/28 16:13:52  nealk
 * Turned off debugging code.
 *
 * Revision 1.1005  1997/02/14 19:18:37  nealk
 * More new/delete commenting.
 *
 * Revision 1.1004  1997/02/14 19:07:30  nealk
 * Added new/delete comments.
 * Played with DPMTA.
 *
 * Revision 1.1003  1997/02/13 16:17:12  ari
 * Intermediate debuging commit - working to fix deep bug in migration?
 *
 * Revision 1.1002  1997/02/11 18:51:41  ari
 * Modified with #ifdef DPMTA to safely eliminate DPMTA codes
 * fixed non-buffering of migration msgs
 * Migration works on multiple processors
 *
 * Revision 1.1001  1997/02/07 22:52:14  jim
 * Eliminated use of nAtomBased and uninitialized memory reads.
 *
 * Revision 1.1000  1997/02/06 15:58:02  ari
 * Resetting CVS to merge branches back into the main trunk.
 * We will stick to main trunk development as suggested by CVS manual.
 * We will set up tags to track fixed points of development/release
 * as suggested by CVS manual - all praise the CVS manual.
 *
 * Revision 1.779  1997/02/06 15:53:00  ari
 * Updating Revision Line, getting rid of branches
 *
 * Revision 1.778.2.1  1997/02/06 02:35:17  jim
 * Implemented periodic boundary conditions - may not work with
 * atom migration yet, but doesn't seem to alter calculation,
 * appears to work correctly when turned on.
 * NamdState chdir's to same directory as config file in argument.
 *
 * Revision 1.778  1997/01/28 00:30:15  ari
 * internal release uplevel to 1.778
 *
 * Revision 1.777  1997/01/17 19:35:48  ari
 * Internal CVS leveling release.  Start development code work
 * at 1.777.1.1.
 *
 * Revision 1.12  1996/12/26 22:26:54  nealk
 * Corrected one of the communications bugs for +p2 -- comm wasn't initialized
 * at the right time.
 *
 * Revision 1.11  1996/12/12 08:57:17  jim
 * added MapDistribMsg packing / unpacking routines
 *
 * Revision 1.10  1996/11/30 01:27:34  jim
 * switched to realistic ComputeType definitions
 *
 * Revision 1.9  1996/11/30 00:32:52  jim
 * added ComputeMgr friend and storage of ComputeType
 *
 * Revision 1.8  1996/11/22 01:02:18  ari
 * *** empty log message ***
 *
 * Revision 1.7  1996/11/07 18:06:56  nealk
 * Changed endl to \n
 *
 * Revision 1.6  1996/11/07 17:55:29  nealk
 * Modified use of DebugM
 *
 * Revision 1.5  1996/11/07 15:21:38  nealk
 * Changed to use DebugM rather than CPrint.
 *
 * Revision 1.4  1996/11/01 21:20:45  ari
 * *** empty log message ***
 *
 * Revision 1.3  1996/10/29 23:35:27  ari
 * *** empty log message ***
 *
 * Revision 1.2  1996/08/16 21:41:11  brunner
 * *** empty log message ***
 *
 * Revision 1.1  1996/08/16 20:43:53  brunner
 * Initial revision
 *
 * Revision 1.1  1996/08/03 20:08:09  brunner
 * Initial revision
 *
 ***************************************************************************/
