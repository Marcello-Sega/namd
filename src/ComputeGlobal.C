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
#include "ComputeGlobalMsgs.h"
#include "PatchMgr.h"
#include "Molecule.h"
#include "ReductionMgr.h"
#include "ComputeMgr.h"
#include "ComputeMgr.decl.h"
#include "SimParameters.h"
#include <stdio.h>

//#define DEBUGM
#define MIN_DEBUG_LEVEL 1
#include "Debug.h"

// CLIENTS

ComputeGlobal::ComputeGlobal(ComputeID c, ComputeMgr *m)
	: ComputeHomePatches(c)
{
  DebugM(3,"Constructing client\n");
    
  ComputeGlobalMaster *master;
  master = 0;
  masterlist.item(0).master = master; 
  comm = m;
}

ComputeGlobal::~ComputeGlobal()
{
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
  DebugM(3,"Receiving configure on client\n");
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
  DebugM(3,"Receiving results (" << msg->aid.size() << " forces, "
	 << msg->newgdef.size() << " new group atoms) on client\n");

  if (msg->tag < 0 || msg->tag >= masterlist.size()) {
    NAMD_die("ComputeGlobal: got ResultsMsg for nonexistent master!");
  }
  MasterConfig &mc = masterlist[msg->tag];

  // set the forces only if we aren't going to resend the data
  int setForces = !msg->resendCoordinates;

  if(setForces) { // we are requested to 
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
      DebugM(1,"processing atom "<<(*a)<<", F="<<(*f2)<<"...\n");
      /* XXX if (*a) is out of bounds here we get a segfault */
      LocalID localID = atomMap->localID(*a);
      if ( localID.pid == notUsed || ! f[localID.pid] ) continue;
      f[localID.pid][localID.index] += (*f2);
    }
    DebugM(1,"done with the loop\n");

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
    DebugM(1,"done with the groups\n");

    for (ap = ap.begin(); ap != ap.end(); ap++) {
      (*ap).forceBox->close(&((*ap).r));
    }

    delete [] f;
  }
  // done setting the forces

  // Get reconfiguration if present
  if ( msg->reconfig ) configure(msg->tag, msg->newaid, msg->newgdef);

  // send another round of data if requested

  if(msg->resendCoordinates) {
    // XXX is this NAMD_die necessary?  Why should it be?
    if(!masterlist[msg->tag].configured)
      NAMD_die("Master requested data before sending reconfig.");

    DebugM(3,"Sending requested data right away\n");
    sendData(msg->tag);
  }

  delete msg;
  DebugM(3,"Done processing results\n");
}

void ComputeGlobal::doWork()
{
  DebugM(2,"doWork\n");
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
  DebugM(2,"done with doWork\n");
}

void ComputeGlobal::sendData(int tag)
{
  DebugM(2,"sendData\n");
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

