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
#include "ComputeMisc.h"
#include "ComputeTcl.h"
//---these include files are needed for ComputeFreeEnergy---
#include <string.h>
#include <strstream.h>
#include "InfoStream.h"
#include "FreeEnergyEnums.h"
#include "FreeEnergyAssert.h"
#include "FreeEnergyGroup.h"
#include "Vector.h"
#include "FreeEnergyVector.h"
#include "FreeEnergyRestrain.h"
#include "FreeEnergyRMgr.h"
#include "FreeEnergyLambda.h"
#include "FreeEnergyLambdMgr.h"
#include "ComputeFreeEnergy.h"
#include "FreeEnergyParse.h"
//----------------------------------------------------------
#include "ComputeIMD.h"
#include "PatchMgr.h"
#include "Molecule.h"
#include "ReductionMgr.h"
#include "ComputeMgr.h"
#include "ComputeMgr.decl.h"
#include "SimParameters.h"
#include <stdio.h>

//#define DEBUGM
#define MIN_DEBUG_LEVEL 4
#include "Debug.h"


// PASS-THROUGH

void ComputeGlobal::recvData(ComputeGlobalDataMsg *msg) {
  if ( master ) {
    master->recvData(msg);
  }
  else NAMD_die("ComputeGlobal::master is NULL!");
}


// CLIENTS

ComputeGlobal::ComputeGlobal(ComputeID c, ComputeMgr *m)
	: ComputeHomePatches(c)
{
  DebugM(3,"Constructing client\n");
  if ( CkMyPe() ) master = 0;
  else {
    SimParameters * simParams = Node::Object()->simParameters;
    if (simParams->IMDon) master = new ComputeIMD(this);
    else if ( simParams->tclForcesOn ) master = new ComputeTcl(this);
    else if ( simParams->miscForcesOn ) master = new ComputeMisc(this);
    else if ( simParams->freeEnergyOn ) master = new ComputeFreeEnergy(this);
    else NAMD_die("Internal error in ComputeGlobal::ComputeGlobal");
  }
  comm = m;
  configured = 0;
}

ComputeGlobal::~ComputeGlobal()
{
  delete master;
}

void ComputeGlobal::configure(AtomIDList newaid, AtomIDList newgdef) {
  DebugM(4,"Receiving configuration (" << newaid.size() <<
	" atoms and " << newgdef.size() << " atoms/groups) on client\n");

  // store data
  aid = newaid;
  gdef = newgdef;

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

  configured = 1;
}

void ComputeGlobal::recvConfig(ComputeGlobalConfigMsg *msg) {
  configure(msg->aid,msg->gdef);
  delete msg;
  sendData();
}

void ComputeGlobal::recvResults(ComputeGlobalResultsMsg *msg) {
  DebugM(3,"Receiving results (" << msg->aid.size() << " forces) on client\n");

  // Store forces to patches
  PatchMap *patchMap = PatchMap::Object();
  int numPatches = patchMap->numPatches();
  AtomMap *atomMap = AtomMap::Object();
  ResizeArrayIter<PatchElem> ap(patchList);
  Force **f = new Force*[numPatches];
  for ( int i = 0; i < numPatches; ++i ) f[i] = 0;

  for (ap = ap.begin(); ap != ap.end(); ap++) {
    (*ap).r = (*ap).forceBox->open();
    f[(*ap).patchID] = (*ap).r->f[Results::normal];
  }

  AtomIDList::iterator a = msg->aid.begin();
  AtomIDList::iterator a_e = msg->aid.end();
  ForceList::iterator f2 = msg->f.begin();
  for ( ; a != a_e; ++a, ++f2 ) {
    LocalID localID = atomMap->localID(*a);
    if ( localID.pid == notUsed || ! f[localID.pid] ) continue;
    f[localID.pid][localID.index] += (*f2);
  }

  // calculate forces for atoms in groups
  Molecule *mol = Node::Object()->molecule;
  AtomIDList::iterator g_i, g_e;
  g_i = gdef.begin(); g_e = gdef.end();
  ResizeArray<BigReal>::iterator gm_i = gmass.begin();
  ForceList::iterator gf_i = msg->gforce.begin();
  for ( ; g_i != g_e; ++g_i, ++gm_i, ++gf_i ) {
    Vector accel = (*gf_i) / (*gm_i);
    for ( ; *g_i != -1; ++g_i ) {
      LocalID localID = atomMap->localID(*g_i);
      if ( localID.pid == notUsed || ! f[localID.pid] ) continue;
      f[localID.pid][localID.index] += accel * mol->atommass(*g_i);
    }
  }

  for (ap = ap.begin(); ap != ap.end(); ap++) {
    (*ap).forceBox->close(&((*ap).r));
  }

  // Get reconfiguration if present
  if ( msg->reconfig ) configure(msg->newaid, msg->newgdef);

  delete msg;
}

void ComputeGlobal::doWork()
{
  if ( configured ) sendData();
  else {
    // send message to check in
    ComputeGlobalDataMsg *msg = new ComputeGlobalDataMsg;
    comm->sendComputeGlobalData(msg);
  }
}

void ComputeGlobal::sendData()
{
  // Get positions from patches
  PatchMap *patchMap = PatchMap::Object();
  int numPatches = patchMap->numPatches();
  AtomMap *atomMap = AtomMap::Object();
  const Lattice & lattice = patchList[0].p->lattice;
  ResizeArrayIter<PatchElem> ap(patchList);
  Position **x = new Position*[numPatches];
  Transform **t = new Transform*[numPatches];
  for ( int i = 0; i < numPatches; ++i ) { x[i] = 0; t[i] = 0; }

  for (ap = ap.begin(); ap != ap.end(); ap++) {
    x[(*ap).patchID] = (*ap).positionBox->open();
    t[(*ap).patchID] = (*ap).p->getTransformList().begin();
    AtomProperties *a = (*ap).atomBox->open();
    (*ap).atomBox->close(&a);
  }

  ComputeGlobalDataMsg *msg = new  ComputeGlobalDataMsg;

  AtomIDList::iterator a = aid.begin();
  AtomIDList::iterator a_e = aid.end();
  for ( ; a != a_e; ++a ) {
    LocalID localID = atomMap->localID(*a);
    if ( localID.pid == notUsed || ! x[localID.pid] ) continue;
    msg->aid.add(*a);
    Position x_orig = x[localID.pid][localID.index];
    Transform trans = t[localID.pid][localID.index];
    msg->p.add(lattice.reverse_transform(x_orig,trans));
  }

  // calculate group centers of mass
  Molecule *mol = Node::Object()->molecule;
  AtomIDList::iterator g_i, g_e;
  g_i = gdef.begin(); g_e = gdef.end();
  ResizeArray<BigReal>::iterator gm_i = gmass.begin();
  for ( ; g_i != g_e; ++g_i, ++gm_i ) {
    Vector com;
    for ( ; *g_i != -1; ++g_i ) {
      LocalID localID = atomMap->localID(*g_i);
      if ( localID.pid == notUsed || ! x[localID.pid] ) continue;
      Position x_orig = x[localID.pid][localID.index];
      Transform trans = t[localID.pid][localID.index];
      com += lattice.reverse_transform(x_orig,trans) * mol->atommass(*g_i);
    }
    com /= *gm_i;
    msg->gcom.add(com);
  }

  for (ap = ap.begin(); ap != ap.end(); ap++) {
    (*ap).positionBox->close(&(x[(*ap).patchID]));
  }
  delete [] x;
  delete [] t;

  DebugM(3,"Sending data (" << msg->aid.size() << " positions) on client\n");

  comm->sendComputeGlobalData(msg);
}

