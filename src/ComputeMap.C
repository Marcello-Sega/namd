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

static char ident[] = "@(#)$Header: /home/cvs/namd/cvsroot/namd2/src/ComputeMap.C,v 1.1003 1997/02/13 16:17:12 ari Exp $";

#include <stdlib.h>
#include <stdio.h>

#include "chare.h"
#include "ckdefs.h"
#include "c++interface.h"

#include "ComputeMap.h"

#define MIN_DEBUG_LEVEL 3
#define DEBUGM 1
#include "Debug.h"

ComputeMap *ComputeMap::_instance = 0;

ComputeMap *ComputeMap::Instance() {
  if (_instance == 0) {
    _instance = new ComputeMap;
  }
  return _instance;
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
	delete[] computeData[i].pids;
    }
    delete[] computeData;
  }    
}


#undef PACK
#define PACK(type,data) { *((type*)b) = data; b += sizeof(type); }

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
  char * const buffer = new char[size];

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
#define UNPACK(type,data) { data = *((type*)b); b += sizeof(type); }

void ComputeMap::unpack (void *in)
{
  DebugM(4,"Unpacking ComputeMap\n");
  int i,j;
  char *b = (char*)in;
  UNPACK(int,nPatchBased);
  UNPACK(int,nAtomBased);
  UNPACK(int,nComputes);
  UNPACK(int,nAllocated);
  computeData = new ComputeData[nAllocated];
  for(i=0;i<nComputes;++i)
  {
    UNPACK(ComputeData,computeData[i]);
    computeData[i].pids = new PatchRec[computeData[i].numPidsAllocated];
    for(j=0;j<computeData[i].numPidsAllocated;++j)
      UNPACK(PatchRec,computeData[i].pids[j]);
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
int ComputeMap::allocateCids(int n)
{
  int i;

  if (computeData != NULL)
  {
    int i;
    for(i=0; i<nComputes; i++)
    {
      if (computeData[i].pids != NULL)
      {
	delete[] computeData[i].pids;
	computeData[i].pids=0;
      }
      delete computeData;

      computeData = NULL;
    }
  }
  nComputes = nPatchBased = nAtomBased = 0;

  computeData = new ComputeData[n];
  nAllocated = n;
  for(i=0; i<n; i++)
  {
    computeData[i].node=0;
    computeData[i].patchBased=false;
    computeData[i].numPids=0;
    computeData[i].pids=NULL;
    computeData[i].compute = NULL;
  }
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
  computeData[cid].pids = new PatchRec[maxPids];
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
  DebugM(1,"---------------------------------------");
  DebugM(1,"---------------------------------------\n");

  DebugM(1,"nComputes = " << nComputes << '\n');
  DebugM(1,"nPatchBased = " << nPatchBased << '\n');
  DebugM(1,"nAtomBased = " << nAtomBased << '\n');
  DebugM(1,"nAllocated = " << nComputes << '\n');
  for(int i=0; i < nComputes; i++)
  {
    DebugM(1,"Compute " << i << '\n');
    DebugM(1,"  node = " << computeData[i].node << '\n');
    DebugM(1,"  patchBased = " << computeData[i].patchBased << '\n');
    DebugM(1,"  numPids = " << computeData[i].numPids << '\n');
    DebugM(1,"  numPidsAllocated = " << computeData[i].numPidsAllocated << '\n');
    for(int j=0; j < computeData[i].numPids; j++)
    {
      DebugM(1,computeData[i].pids[j].pid);
      if (!((j+1) % 6))
	DebugM(1,'\n');
    }
    DebugM(1,"\n---------------------------------------");
    DebugM(1,"---------------------------------------\n");

  }
}
/***************************************************************************
 * RCS INFORMATION:
 *
 *	$RCSfile: ComputeMap.C,v $
 *	$Author: ari $	$Locker:  $		$State: Exp $
 *	$Revision: 1.1003 $	$Date: 1997/02/13 16:17:12 $
 *
 ***************************************************************************
 * REVISION HISTORY:
 *
 * $Log: ComputeMap.C,v $
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
