/**
***  Copyright (c) 1995, 1996, 1997, 1998, 1999, 2000 by
***  The Board of Trustees of the University of Illinois.
***  All rights reserved.
**/

/*
   Forwards atoms to master node for force evaluation.
*/

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
#if !defined(WIN32) || defined(__CYGWIN__)
#include <strstream.h>
#else
#include <strstrea.h>
#endif
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
#include "ComputeSMD.h"
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
  int tag = msg->tag;
  if (tag >= 0 && tag < masterlist.size() && masterlist[tag].master != 0) {
    masterlist[tag].master->recvData(msg);
  } else {
    char buf[1024];
    sprintf(buf, "ComputeGlobal::received DataMsg for nonexistent master!  tag=%d, size=%d", msg->tag, masterlist.size());
    NAMD_die(buf);
  }
}
    
// CLIENTS

ComputeGlobal::ComputeGlobal(ComputeID c, ComputeMgr *m)
	: ComputeHomePatches(c)
{
  DebugM(3,"Constructing client\n");
  ComputeGlobalMaster *master;
  if ( !CkMyPe() ) { 
    SimParameters * simParams = Node::Object()->simParameters;
    if (simParams->IMDon) master = new ComputeIMD(m);
    else if ( simParams->tclForcesOn ) master = new ComputeTcl(m);
    else if ( simParams->miscForcesOn ) master = new ComputeMisc(m);
    else if ( simParams->freeEnergyOn ) master = new ComputeFreeEnergy(m);
    else if ( simParams->SMDOn ) master = new ComputeSMD(m);
    else NAMD_die("Internal error in ComputeGlobal::ComputeGlobal");
 
    master->set_tag(0); // Hard-coded until we can have multiple masters
  } else {
    master = 0;
  }
  masterlist.item(0).master = master; 
  comm = m;
}

ComputeGlobal::~ComputeGlobal()
{
  for (int i=0; i<masterlist.size(); i++) 
    delete masterlist[i].master;
}

void ComputeGlobal::configure(int tag, AtomIDList newaid, AtomIDList newgdef) {
  DebugM(4,"Receiving configuration (" << newaid.size() <<
	" atoms and " << newgdef.size() << " atoms/groups) on client\n");

  MasterConfig &mc = masterlist[tag]; 

  // store data
  mc.aid = newaid;
  mc.gdef = newgdef;

  // calculate group masses
  Molecule *mol = Node::Object()->molecule;
  mc.gmass.resize(0);
  AtomIDList::iterator g_i, g_e;
  g_i = mc.gdef.begin(); g_e = mc.gdef.end();
  for ( ; g_i != g_e; ++g_i ) {
    BigReal mass = 0;
    for ( ; *g_i != -1; ++g_i ) {
      mass += mol->atommass(*g_i);
    }
    mc.gmass.add(mass);
  }

  mc.configured = 1;
}

void ComputeGlobal::recvConfig(ComputeGlobalConfigMsg *msg) {
  int tag = msg->tag;
  if (tag >= 0 && tag < masterlist.size()) { 
    configure(tag, msg->aid,msg->gdef);
    delete msg;
    sendData(tag);
  } else {
    NAMD_die("ComputeGlobal: received ConfigMsg for nonexistent master!");
  }
}

void ComputeGlobal::recvResults(ComputeGlobalResultsMsg *msg) {
  DebugM(3,"Receiving results (" << msg->aid.size() << " forces) on client\n");

  if (msg->tag < 0 || msg->tag >= masterlist.size()) {
    NAMD_die("ComputeGlobal: got ResultsMsg for nonexistent master!");
  }
  MasterConfig &mc = masterlist[msg->tag];
  
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
  int q=0;
  for ( ; a != a_e; ++a, ++f2 ) {
    LocalID localID = atomMap->localID(*a);
    if ( localID.pid == notUsed || ! f[localID.pid] ) continue;
    f[localID.pid][localID.index] += (*f2);
  }

  // calculate forces for atoms in groups
  Molecule *mol = Node::Object()->molecule;
  AtomIDList::iterator g_i, g_e;
  g_i = mc.gdef.begin(); g_e = mc.gdef.end();
  ResizeArray<BigReal>::iterator gm_i = mc.gmass.begin();
  ForceList::iterator gf_i = msg->gforce.begin();
  //iout << iDEBUG << "recvResults\n" << endi;
  for ( ; g_i != g_e; ++g_i, ++gm_i, ++gf_i ) {
    //iout << iDEBUG << *gf_i << '\n' << endi;
    Vector accel = (*gf_i) / (*gm_i);
    for ( ; *g_i != -1; ++g_i ) {
      //iout << iDEBUG << *g_i << '\n' << endi;
      LocalID localID = atomMap->localID(*g_i);
      if ( localID.pid == notUsed || ! f[localID.pid] ) continue;
      f[localID.pid][localID.index] += accel * mol->atommass(*g_i);
    }
  }

  for (ap = ap.begin(); ap != ap.end(); ap++) {
    (*ap).forceBox->close(&((*ap).r));
  }

  // Get reconfiguration if present
  if ( msg->reconfig ) configure(msg->tag, msg->newaid, msg->newgdef);

  delete [] f;
  delete msg;
}

void ComputeGlobal::doWork()
{
  for (int i=0; i<masterlist.size(); i++) {
    if (masterlist[i].configured) {
      sendData(i);
    } else {
      // send message to check in
      ComputeGlobalDataMsg *msg = new ComputeGlobalDataMsg;
      msg->tag = i;
      comm->sendComputeGlobalData(msg);
    }
  }
}

void ComputeGlobal::sendData(int tag)
{
  // Get positions from patches
  PatchMap *patchMap = PatchMap::Object();
  int numPatches = patchMap->numPatches();
  AtomMap *atomMap = AtomMap::Object();
  const Lattice & lattice = patchList[0].p->lattice;
  ResizeArrayIter<PatchElem> ap(patchList);
  CompAtom **x = new CompAtom*[numPatches];
  FullAtom **t = new FullAtom*[numPatches];
  for ( int i = 0; i < numPatches; ++i ) { x[i] = 0; t[i] = 0; }

  for (ap = ap.begin(); ap != ap.end(); ap++) {
    x[(*ap).patchID] = (*ap).positionBox->open();
    t[(*ap).patchID] = (*ap).p->getAtomList().begin();
  }

  ComputeGlobalDataMsg *msg = new  ComputeGlobalDataMsg;
  msg->tag = tag;
  MasterConfig &mc = masterlist[tag];

  AtomIDList::iterator a = mc.aid.begin();
  AtomIDList::iterator a_e = mc.aid.end();
  for ( ; a != a_e; ++a ) {
    LocalID localID = atomMap->localID(*a);
    if ( localID.pid == notUsed || ! x[localID.pid] ) continue;
    msg->aid.add(*a);
    Position x_orig = x[localID.pid][localID.index].position;
    Transform trans = t[localID.pid][localID.index].transform;
    msg->p.add(lattice.reverse_transform(x_orig,trans));
  }

  // calculate group centers of mass
  Molecule *mol = Node::Object()->molecule;
  AtomIDList::iterator g_i, g_e;
  g_i = mc.gdef.begin(); g_e = mc.gdef.end();
  ResizeArray<BigReal>::iterator gm_i = mc.gmass.begin();
  for ( ; g_i != g_e; ++g_i, ++gm_i ) {
    Vector com(0,0,0);
    for ( ; *g_i != -1; ++g_i ) {
      LocalID localID = atomMap->localID(*g_i);
      if ( localID.pid == notUsed || ! x[localID.pid] ) continue;
      Position x_orig = x[localID.pid][localID.index].position;
      Transform trans = t[localID.pid][localID.index].transform;
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

