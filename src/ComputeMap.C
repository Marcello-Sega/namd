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

static char ident[] = "@(#)$Header: /home/cvs/namd/cvsroot/namd2/src/ComputeMap.C,v 1.5 1996/11/07 15:21:38 nealk Exp $";

#include <stdlib.h>
#include <stdio.h>

#include "chare.h"
#include "ckdefs.h"
#include "c++interface.h"

#include "ComputeMap.h"

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
    int i;
    for(i=0; i<nComputes; i++)
    {
      if (computeData[i].pids != NULL)
      {
	delete computeData[i].pids;
	computeData[i].pids=0;
      }
      delete computeData;

      computeData = NULL;
    }
  }    
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
    return computeData[cid].pids[i];
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
	delete computeData[i].pids;
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

  if (angleForceType == type)
  {
    computeData[cid].patchBased = false;
    nAtomBased++;
  }
  else
  {
    computeData[cid].patchBased = true;
    nPatchBased++;
  }

  computeData[cid].numPids = 0;
  computeData[cid].pids = new int[maxPids];
  computeData[cid].numPidsAllocated = maxPids;

  return cid;
}

//----------------------------------------------------------------------
int ComputeMap::newPid(ComputeID cid, PatchID pid)
{
  if (!computeData)
    return -1;                   // Have to allocate first

  if ((cid < 0) || (cid >= nComputes))
    return -1;                   // Have to store the cid first

  if (computeData[cid].numPids == computeData[cid].numPidsAllocated)
    return -1;                   // Out of space for dependents

  computeData[cid].pids[computeData[cid].numPids]=pid;
  computeData[cid].numPids++;
  return 0;
}

//----------------------------------------------------------------------
void ComputeMap::printComputeMap(void)
{
  DebugM(1,"---------------------------------------");
  DebugM(1,"---------------------------------------\n");

  DebugM(1,"nComputes = %d\n",nComputes);
  DebugM(1,"nPatchBased = %d\n",nPatchBased);
  DebugM(1,"nAtomBased = %d\n",nAtomBased);
  DebugM(1,"nAllocated = %d\n",nComputes);
  for(int i=0; i < nComputes; i++)
  {
    DebugM(1,"Compute %d\n",i);
    DebugM(1,"  node = %d\n",computeData[i].node);
    DebugM(1,"  patchBased = %d\n",computeData[i].patchBased);
    DebugM(1,"  numPids = %d\n",computeData[i].numPids);
    DebugM(1,"  numPidsAllocated = %d\n",computeData[i].numPidsAllocated);
    for(int j=0; j < computeData[i].numPids; j++)
    {
      DebugM(1," %10d ",computeData[i].pids[j]);
      if (!((j+1) % 6))
	DebugM(1,"\n");
    }
    DebugM(1,"\n---------------------------------------");
    DebugM(1,"---------------------------------------\n");

  }
}
/***************************************************************************
 * RCS INFORMATION:
 *
 *	$RCSfile: ComputeMap.C,v $
 *	$Author: nealk $	$Locker:  $		$State: Exp $
 *	$Revision: 1.5 $	$Date: 1996/11/07 15:21:38 $
 *
 ***************************************************************************
 * REVISION HISTORY:
 *
 * $Log: ComputeMap.C,v $
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
