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

#include "Node.h"
#include "ComputeNonbondedUtil.h"
#include "ComputeNonbondedSelf.h"
#include "ComputeNonbondedPair.h"
#include "ComputeAngles.h"
#include "ComputeDihedrals.h"
#include "ComputeImpropers.h"
#include "ComputeBonds.h"
#include "ComputeNonbondedExcl.h"
#define MIN_DEBUG_LEVEL 4
#define DEBUGM
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
  DebugM(2,"createComputes 0\n");
  Node *node = Node::Object();
  int myNode = node->myid();

  int numNonbondedSelf = 0;
  int numNonbondedPair = 0;

  DebugM(1,"---------------------------------------\n");
  DebugM(1,"---------------------------------------\n");

  DebugM(3,"nComputes = " << map->nComputes << '\n');
  DebugM(3,"nPatchBased = " << map->nPatchBased << '\n');
  DebugM(3,"nAtomBased = " << map->nAtomBased << '\n');
  DebugM(3,"nAllocated = " << map->nComputes << '\n');
  DebugM(2,"createComputes 1: looping " << map->nComputes << "\n");
  for(int i=0; i < map->nComputes; i++)
  {
    if ( map->computeData[i].node != myNode ) continue;
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

  DebugM(2,"createComputes 2: looping " << i << "on type: " << map->computeData[i].type << "\n");
    switch ( map->computeData[i].type )
    {
      case computeNonbondedSelfType:
	c = new ComputeNonbondedSelf(i,map->computeData[i].pids[0]);
	DebugM(3,"ComputeNonbondedSelf created.\n");
	++numNonbondedSelf;
	map->registerCompute(i,c);
	c->mapReady();
	break;
      case computeNonbondedPairType:
	pid2[0] = map->computeData[i].pids[0];
	pid2[1] = map->computeData[i].pids[1];
	c = new ComputeNonbondedPair(i,pid2);
	DebugM(3,"ComputeNonbondedPair created.\n");
	++numNonbondedPair;
	map->registerCompute(i,c);
	c->mapReady();
	break;
      case computeNonbondedExclType:
	c = new ComputeNonbondedExcls(i);
	DebugM(3,"ComputeNonbondedExcls created.\n");
	map->registerCompute(i,c);
	DebugM(3,"ComputeNonbondedExcls registered.\n");
	c->mapReady();
	DebugM(3,"ComputeNonbondedExcls ready.\n");
	break;
      case computeBondsType:
	c = new ComputeBonds(i);
	DebugM(3,"ComputeBonds created.\n");
	map->registerCompute(i,c);
	c->mapReady();
	break;
      case computeAnglesType:
	c = new ComputeAngles(i);
	DebugM(3,"ComputeAngles created.\n");
	map->registerCompute(i,c);
	c->mapReady();
	break;
      case computeDihedralsType:
	c = new ComputeDihedrals(i);
	DebugM(3,"ComputeDihedrals created.\n");
	map->registerCompute(i,c);
	c->mapReady();
	break;
      case computeImpropersType:
	c = new ComputeImpropers(i);
	DebugM(3,"ComputeImpropers created.\n");
	map->registerCompute(i,c);
	c->mapReady();
	break;
      default:
	DebugM(10,"Unknown compute type not created!\n");
    }

  }

  DebugM(2,"createComputes 5: done looping\n");
  DebugM(4, numNonbondedSelf << " ComputeNonbondedSelf created\n");
  DebugM(4, numNonbondedPair << " ComputeNonbondedPair created\n");

  // This doesn't really have to be here, but the output makes more sense.
  // It does have to happen after the molecule has been created.
  ComputeNonbondedUtil::select();

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
 *	$Author: ari $	$Locker:  $		$State: Exp $
 *	$Revision: 1.777 $	$Date: 1997/01/17 19:35:49 $
 *
 ***************************************************************************
 * REVISION HISTORY:
 *
 *
 ***************************************************************************/
