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
#include "PatchMap.h"
#include "PatchMap.inl"
#include "AtomMap.h"
#include "ComputeGlobal.h"
#include "ComputeGlobalMaster.h"
#include "ComputeGlobalMsgs.h"
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
#include "ComputeMDComm.h"
#include "PatchMgr.h"
#include "Molecule.h"
#include "ReductionMgr.h"
#include "ComputeMgr.h"
#include "ComputeMgr.top.h"
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
  if ( CMyPe() ) master = 0;
  else {
    SimParameters * simParams = Node::Object()->simParameters;
    if ( simParams->tclForcesOn ) master = new ComputeTcl(this);
    else if ( simParams->freeEnergyOn ) master = new ComputeFreeEnergy(this);
#ifdef MDCOMM
    else if ( simParams->vmdFrequency != -1 ) master = new ComputeMDComm(this);
#endif
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
    ComputeGlobalDataMsg *msg =
	new (MsgIndex(ComputeGlobalDataMsg)) ComputeGlobalDataMsg;
    comm->sendComputeGlobalData(msg);
  }
}

void ComputeGlobal::sendData()
{
  // Get positions from patches
  PatchMap *patchMap = PatchMap::Object();
  int numPatches = patchMap->numPatches();
  AtomMap *atomMap = AtomMap::Object();
  ResizeArrayIter<PatchElem> ap(patchList);
  Position **x = new Position*[numPatches];
  for ( int i = 0; i < numPatches; ++i ) x[i] = 0;

  for (ap = ap.begin(); ap != ap.end(); ap++) {
    x[(*ap).patchID] = (*ap).positionBox->open();
    AtomProperties *a = (*ap).atomBox->open();
    (*ap).atomBox->close(&a);
  }

  ComputeGlobalDataMsg *msg =
	new (MsgIndex(ComputeGlobalDataMsg)) ComputeGlobalDataMsg;

  AtomIDList::iterator a = aid.begin();
  AtomIDList::iterator a_e = aid.end();
  for ( ; a != a_e; ++a ) {
    LocalID localID = atomMap->localID(*a);
    if ( localID.pid == notUsed || ! x[localID.pid] ) continue;
    msg->aid.add(*a);
    msg->p.add(x[localID.pid][localID.index]);
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
      com += x[localID.pid][localID.index] * mol->atommass(*g_i);
    }
    com /= *gm_i;
    msg->gcom.add(com);
  }

  for (ap = ap.begin(); ap != ap.end(); ap++) {
    (*ap).positionBox->close(&(x[(*ap).patchID]));
  }
  delete [] x;

  DebugM(3,"Sending data (" << msg->aid.size() << " positions) on client\n");

  comm->sendComputeGlobalData(msg);
}


/***************************************************************************
 * RCS INFORMATION:
 *
 *	$RCSfile $
 *	$Author $	$Locker:  $		$State: Exp $
 *	$Revision: 1.11 $	$Date: 1998/05/21 22:43:23 $
 *
 ***************************************************************************
 * REVISION HISTORY:
 *
 * $Log: ComputeGlobal.C,v $
 * Revision 1.11  1998/05/21 22:43:23  hurwitz
 * initial check in of fixed and forcing restraints
 * -Dave Hurwitz
 *
 * Revision 1.10  1998/04/30 04:53:22  jim
 * Added forces from MDComm and other improvements to ComputeGlobal.
 *
 * Revision 1.9  1998/02/16 00:47:27  jim
 * Added code to actually distribute group forces to atoms.
 *
 * Revision 1.8  1998/02/16 00:23:18  jim
 * Added atom group centers of mass to Tcl interface.
 *
 * Revision 1.7  1998/02/10 06:45:09  jim
 * Added class ComputeFreeEnergy.
 *
 * Revision 1.6  1998/02/10 05:35:02  jim
 * Split ComputeGlobal into different classes and files.
 * Switched globalForces and globalForcesTcl to tclForces and tclForcesScript.
 * Added (soon to be used) freeEnergy and freeEnergyConfig.
 *
 * Revision 1.5  1998/01/15 04:58:45  jim
 * Corrected "friend foo" to "friend class foo".
 *
 * Revision 1.4  1998/01/10 23:57:25  jim
 * Added loadmasses command to TCL.
 *
 * Revision 1.3  1998/01/06 05:41:26  jim
 * Added tclx library.
 *
 * Revision 1.2  1997/12/26 23:10:43  milind
 * Made namd2 to compile, link and run under linux. Merged Templates and src
 * directoriies, and removed separate definition and declaration files for
 * templates.
 *
 * Revision 1.1  1997/12/19 23:48:45  jim
 * Added Tcl interface for calculating forces.
 *
 *
 ***************************************************************************/
