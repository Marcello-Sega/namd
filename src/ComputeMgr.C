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

#define MIN_DEBUG_LEVEL 4
//#define DEBUGM
#include "Debug.h"

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
#include "ComputeFullDirect.h"
#include "ComputeDPMTA.h"
#include "ComputeCylindricalBC.h"

ComputeMgr::ComputeMgr(InitMsg *msg)
{
  delete msg;
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
    if ( ! ( i % 100 ) )
    {
      DebugM(4,"Created " << i << " compute objects so far.\n");
    }
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
    int trans2[2];

  DebugM(2,"createComputes 2: looping " << i << "on type: " << map->computeData[i].type << "\n");
    switch ( map->computeData[i].type )
    {
      case computeNonbondedSelfType:
	c = new ComputeNonbondedSelf(i,map->computeData[i].pids[0].pid); // unknown delete
	DebugM(3,"ComputeNonbondedSelf created.\n");
	++numNonbondedSelf;
	map->registerCompute(i,c);
	c->initialize();
	break;
      case computeNonbondedPairType:
	pid2[0] = map->computeData[i].pids[0].pid;
	trans2[0] = map->computeData[i].pids[0].trans;
	pid2[1] = map->computeData[i].pids[1].pid;
	trans2[1] = map->computeData[i].pids[1].trans;
	c = new ComputeNonbondedPair(i,pid2,trans2); // unknown delete
	DebugM(3,"ComputeNonbondedPair created.\n");
	++numNonbondedPair;
	map->registerCompute(i,c);
	c->initialize();
	break;
      case computeNonbondedExclType:
	c = new ComputeNonbondedExcls(i); // unknown delete
	DebugM(4,"ComputeNonbondedExcls created.\n");
	map->registerCompute(i,c);
	DebugM(3,"ComputeNonbondedExcls registered.\n");
	c->initialize();
	DebugM(3,"ComputeNonbondedExcls ready.\n");
	break;
      case computeBondsType:
	c = new ComputeBonds(i); // unknown delete
	DebugM(4,"ComputeBonds created.\n");
	map->registerCompute(i,c);
	c->initialize();
	break;
      case computeAnglesType:
	c = new ComputeAngles(i); // unknown delete
	DebugM(4,"ComputeAngles created.\n");
	map->registerCompute(i,c);
	c->initialize();
	break;
      case computeDihedralsType:
	c = new ComputeDihedrals(i); // unknown delete
	DebugM(4,"ComputeDihedrals created.\n");
	map->registerCompute(i,c);
	c->initialize();
	break;
      case computeImpropersType:
	c = new ComputeImpropers(i); // unknown delete
	DebugM(4,"ComputeImpropers created.\n");
	map->registerCompute(i,c);
	c->initialize();
	break;
#ifdef DPMTA
      case computeDPMTAType:
	c = new ComputeDPMTA(i); // unknown delete
	DebugM(4,"ComputeDPMTA created.\n");
	map->registerCompute(i,c);
	c->initialize();
	break;
      case computeDPMEType:
	// c = new ComputeDPME(i); // unknown delete
	DebugM(4,"ComputeDPME *NOT* created.\n");
	// map->registerCompute(i,c);
	// c->initialize();
	break;
#endif
      case computeFullDirectType:
	c = new ComputeFullDirect(i); // unknown delete
	DebugM(4,"ComputeFullDirect created.\n");
	map->registerCompute(i,c);
	c->initialize();
	break;
      case computeCylindricalBCType:
	c = new ComputeCylindricalBC(i,map->computeData[i].pids[0].pid); // unknown delete
	DebugM(4,"ComputeCylindricalBC created.\n");
	map->registerCompute(i,c);
	c->initialize();
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


#include "ComputeMgr.bot.h"

/***************************************************************************
 * RCS INFORMATION:
 *
 *	$RCSfile: ComputeMgr.C,v $
 *	$Author: jim $	$Locker:  $		$State: Exp $
 *	$Revision: 1.1006 $	$Date: 1997/03/15 22:15:25 $
 *
 ***************************************************************************
 * REVISION HISTORY:
 *
 *
 ***************************************************************************/
