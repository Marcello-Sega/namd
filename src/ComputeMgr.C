/***************************************************************************/
/*                                                                         */
/*              (C) Copyright 1996 The Board of Trustees of the            */
/*                          University of Illinois                         */
/*                           All Rights Reserved                           */
/*									   */
/***************************************************************************/

/***************************************************************************
 * DESCRIPTION:
 *
 ***************************************************************************/

#include "ckdefs.h"
#include "chare.h"
#include "c++interface.h"

#include "BOCgroup.h"
#include "ComputeMgr.top.h"
#include "ComputeMgr.h"

#include "ComputeNonbondedSelf.h"
#include "ComputeNonbondedPair.h"
#include "ComputeAngles.h"

#define MIN_DEBUG_LEVEL 4
#define DEBUGM 1
#include "Debug.h"

ComputeMgr::ComputeMgr(InitMsg *msg)
{
  ;
}

ComputeMgr::~ComputeMgr(void)
{
  ;
}

void ComputeMgr:: createComputes(ComputeMap *map)
{
  int numNonbondedSelf = 0;
  int numNonbondedPair = 0;

  DebugM(1,"---------------------------------------");
  DebugM(1,"---------------------------------------\n");

  DebugM(1,"nComputes = " << map->nComputes << '\n');
  DebugM(1,"nPatchBased = " << map->nPatchBased << '\n');
  DebugM(1,"nAtomBased = " << map->nAtomBased << '\n');
  DebugM(1,"nAllocated = " << map->nComputes << '\n');
  for(int i=0; i < map->nComputes; i++)
  {
    DebugM(1,"Compute " << i << '\n');
    DebugM(1,"  node = " << map->computeData[i].node << '\n');
    DebugM(1,"  type = " << map->computeData[i].type << '\n');
    DebugM(1,"  patchBased = " << map->computeData[i].patchBased << '\n');
    DebugM(1,"  numPids = " << map->computeData[i].numPids << '\n');
    DebugM(1,"  numPidsAllocated = " << map->computeData[i].numPidsAllocated << '\n');
    for(int j=0; j < map->computeData[i].numPids; j++)
    {
      DebugM(1,"  pid " << map->computeData[i].pids[j] << '\n');
      if (!((j+1) % 6))
	DebugM(1,'\n');
    }
    DebugM(1,"\n---------------------------------------");
    DebugM(1,"---------------------------------------\n");

    Compute *c;
    PatchID pid2[2];

    switch ( map->computeData[i].type )
    {
      case computeNonbondedSelfType:
	c = new ComputeNonbondedSelf(i,map->computeData[i].pids[0]);
	DebugM(3,"ComputeNonbondedSelf created\n");
	++numNonbondedSelf;
	map->registerCompute(i,c);
	c->mapReady();
	break;
      case computeNonbondedPairType:
	pid2[0] = map->computeData[i].pids[0];
	pid2[1] = map->computeData[i].pids[1];
	c = new ComputeNonbondedPair(i,pid2);
	DebugM(3,"ComputeNonbondedPair created\n");
	++numNonbondedPair;
	map->registerCompute(i,c);
	c->mapReady();
	break;
      case computeNonbondedExclType:
	DebugM(2,"ComputeNonbondedExcl would have been created.\n");
	break;
      case computeBondsType:
	DebugM(2,"ComputeBonds would have been created.\n");
	break;
      case computeAnglesType:
	c = new ComputeAngles(i);
	DebugM(3,"ComputeAngles created\n");
	map->registerCompute(i,c);
	c->mapReady();
	break;
      case computeDihedralsType:
	DebugM(2,"ComputeDihedrals would have been created.\n");
	break;
      case computeImpropersType:
	DebugM(2,"ComputeImpropers would have been created.\n");
	break;
      default:
	DebugM(10,"Unknown compute type not created!\n");
    }

  }

  DebugM(4, numNonbondedSelf << " ComputeNonbondedSelf created\n");
  DebugM(4, numNonbondedPair << " ComputeNonbondedPair created\n");

}

void ComputeMgr:: enqueueWork(Compute *compute)
{
  ;
}

#include "ComputeMgr.bot.h"

/***************************************************************************
 * RCS INFORMATION:
 *
 *	$RCSfile: ComputeMgr.C,v $
 *	$Author: jim $	$Locker:  $		$State: Exp $
 *	$Revision: 1.3 $	$Date: 1996/11/30 02:02:11 $
 *
 ***************************************************************************
 * REVISION HISTORY:
 *
 *
 ***************************************************************************/
