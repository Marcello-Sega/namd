/**
***  Copyright (c) 1995, 1996, 1997, 1998, 1999, 2000 by
***  The Board of Trustees of the University of Illinois.
***  All rights reserved.
**/

/*
   Forwards atoms to master node for force evaluation.
*/

#include "Namd.h"
#include "Node.h"
#include "PatchMap.h"
#include "PatchMap.inl"
#include "AtomMap.h"
#include "ComputeGlobal.h"
#include "ComputeGlobalMaster.h"
#include "ComputeGlobalMsgs.h"
#include "PatchMgr.h"
#include "Molecule.h"
#include "ReductionMgr.h"
#include "ComputeMgr.h"
#include "ComputeMgr.decl.h"
#include <stdio.h>

// #define DEBUGM
#define MIN_DEBUG_LEVEL 4
#include "Debug.h"


ComputeGlobalMaster::ComputeGlobalMaster(ComputeGlobal *h) {
  DebugM(3,"Constructing master\n");
  host = h;
  initialized = 0;
  msgcount = 0;
  numWorkingPes = CkNumPes();
  {
    int npatches=(PatchMap::Object())->numPatches();
    if ( numWorkingPes > npatches ) numWorkingPes = npatches;
  }
}

ComputeGlobalMaster::~ComputeGlobalMaster() {
  DebugM(3,"Destructing master\n");
}

void ComputeGlobalMaster::recvData(ComputeGlobalDataMsg *msg) {
  DebugM(3,"Receiving data on master\n");
  // Check initialization and number of messages received
  if ( ! initialized )
  {
    delete msg;
    if ( ++msgcount == numWorkingPes ) {
       msgcount = 0;
       initialize();
       initialized = 1;
    }
  }
  else {
    storedata(msg);
    if ( ++msgcount == numWorkingPes ) {
       msgcount = 0;
       calculate();
       cleardata();
    }
  }
}

void ComputeGlobalMaster::initialize() {
  DebugM(4,"Initializing master\n");

  ComputeGlobalConfigMsg *msg = new ComputeGlobalConfigMsg;

  // Build config here

  // Send config to clients
  host->comm->sendComputeGlobalConfig(msg);
}

void ComputeGlobalMaster::storedata(ComputeGlobalDataMsg *msg) {
  DebugM(3,"Storing data (" << msg->aid.size() << " positions) on master\n");
  AtomIDList::iterator a_i = msg->aid.begin();
  AtomIDList::iterator a_e = msg->aid.end();
  PositionList::iterator p_i = msg->p.begin();
  for ( ; a_i != a_e; ++a_i, ++p_i ) {
    aid.add(*a_i);
    p.add(*p_i);
  }
  if ( ! gcom.size() ) gcom.resize(msg->gcom.size());
  PositionList::iterator c_i, c_e;
  c_i = msg->gcom.begin(); c_e = msg->gcom.end();
  p_i = gcom.begin();
  for ( ; c_i != c_e; ++c_i, ++p_i ) {
    *p_i += *c_i;
  }
  delete msg;
}

void ComputeGlobalMaster::cleardata() {
  aid.resize(0);
  p.resize(0);
  gcom.resize(0);
}

void ComputeGlobalMaster::storedefs(AtomIDList newgdef) {
  // store group definitions
  gdef = newgdef;  // aliases, but original should be deleted

  // calculate group masses
  Molecule *mol = Node::Object()->molecule;
  gmass.resize(0);
  AtomIDList::iterator g_i, g_e;
  g_i = gdef.begin(); g_e = gdef.end();
  for ( ; g_i != g_e; ++g_i ) {
    BigReal mass = 0;
    for ( ; *g_i != -1; ++g_i ) {
      mass += mol->atommass(*g_i);
    }
    gmass.add(mass);
  }
}

void ComputeGlobalMaster::calculate() {
  DebugM(4,"Calculating forces on master\n");

  ComputeGlobalResultsMsg *msg = new ComputeGlobalResultsMsg;
  msg->gforce.resize(gmass.size());

  // Build results here

  // Send results to clients
  DebugM(3,"Sending results (" << msg->aid.size() << " forces) on master\n");
  if ( msg->reconfig ) {
    DebugM(4,"Sending new configuration (" <<
			msg->newaid.size() << " atoms) on master\n");
  }
  host->comm->sendComputeGlobalResults(msg);
}

