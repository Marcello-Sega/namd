/***************************************************************************/
/*       (C) Copyright 1996,1997 The Board of Trustees of the              */
/*                          University of Illinois                         */
/*                           All Rights Reserved                           */
/***************************************************************************/
/***************************************************************************
 * DESCRIPTION:  Forwards atoms to master node for force evaluation.
 *
 ***************************************************************************/

#include "Namd.h"
#include "Node.h"
#include "ComputeFreeEnergy.h"
#include "ComputeGlobal.h"
#include "ComputeGlobalMsgs.h"
#include "Molecule.h"
#include "ReductionMgr.h"
#include "ComputeMgr.h"
#include "ComputeMgr.top.h"
#include <stdio.h>

// #define DEBUGM
#define MIN_DEBUG_LEVEL 4
#include "Debug.h"


ComputeFreeEnergy::ComputeFreeEnergy(ComputeGlobal *h) : ComputeGlobalMaster(h) {
  DebugM(3,"Constructing ComputeFreeEnergy\n");
}

ComputeFreeEnergy::~ComputeFreeEnergy() {
  DebugM(3,"Destructing ComputeFreeEnergy\n");
}


void ComputeFreeEnergy::initialize() {
  DebugM(4,"Initializing master\n");

  ComputeGlobalConfigMsg *msg =
	new (MsgIndex(ComputeGlobalConfigMsg)) ComputeGlobalConfigMsg;

  // Get the path for our script
  char *filename = Node::Object()->configList->find("freeEnergyConfig")->data;

  iout << iDEBUG << "Free energy perturbation - initialize()\n" << endi; 

  // Send config to clients
  host->comm->sendComputeGlobalConfig(msg);
}


void ComputeFreeEnergy::calculate() {
  DebugM(4,"Calculating forces on master\n");

  ComputeGlobalResultsMsg *msg =
	new (MsgIndex(ComputeGlobalResultsMsg)) ComputeGlobalResultsMsg;

  iout << iDEBUG << "Free energy perturbation - calculate()\n" << endi; 

  // Send results to clients
  DebugM(3,"Sending results (" << msg->aid.size() << " forces) on master\n");
  if ( msg->reconfig ) {
    DebugM(4,"Sending new configuration (" <<
			msg->newaid.size() << " atoms) on master\n");
  }
  host->comm->sendComputeGlobalResults(msg);
}


/***************************************************************************
 * RCS INFORMATION:
 *
 *	$RCSfile $
 *	$Author $	$Locker:  $		$State: Exp $
 *	$Revision: 1.1 $	$Date: 1998/02/10 06:45:08 $
 *
 ***************************************************************************
 * REVISION HISTORY:
 *
 * $Log: ComputeFreeEnergy.C,v $
 * Revision 1.1  1998/02/10 06:45:08  jim
 * Added class ComputeFreeEnergy.
 *
 *
 *
 *
 ***************************************************************************/
