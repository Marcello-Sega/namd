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
#include "ComputeGlobalMsgs.h"
#include "PatchMgr.h"
#include "Molecule.h"
#include "ReductionMgr.h"
#include "ComputeMgr.h"
#include "ComputeMgr.top.h"
#include <stdio.h>

#ifdef NAMD_TCL
#include <tcl.h>
#include "TclCommands.h"
#endif

// #define DEBUGM
#define MIN_DEBUG_LEVEL 4
#include "Debug.h"

// MASTER

class ComputeGlobalMaster {
private:
  friend ComputeGlobal;
  ComputeGlobal *host;
  ComputeGlobalMaster(ComputeGlobal *);
  ~ComputeGlobalMaster();
  void recvData(ComputeGlobalDataMsg *);
  int msgcount;
  void initialize();
  int initialized;
  void storedata(ComputeGlobalDataMsg *);
  void cleardata();
  AtomIDList aid;
  PositionList p;
  void calculate();
#ifdef NAMD_TCL
  Tcl_Interp *interp;
  static int Tcl_print(ClientData, Tcl_Interp *, int, char **);
  static int Tcl_addatom(ClientData, Tcl_Interp *, int, char **);
  static int Tcl_reconfig(ClientData, Tcl_Interp *, int, char **);
  static int Tcl_loadcoords(ClientData, Tcl_Interp *, int, char **);
  static int Tcl_addforce(ClientData, Tcl_Interp *, int, char **);
#endif
};


int ComputeGlobalMaster::Tcl_print(ClientData,
	Tcl_Interp *, int argc, char *argv[]) {
  char *msg = Tcl_Merge(argc-1,argv+1);
  CPrintf("TCL: %s\n",msg);
  free(msg);
  return TCL_OK;
}


int ComputeGlobalMaster::Tcl_addatom(ClientData clientData,
	Tcl_Interp *interp, int argc, char *argv[]) {
  if (argc != 2) {
    interp->result = "wrong # args";
    return TCL_ERROR;
  }
  int atomid;
  if (Tcl_GetInt(interp,argv[1],&atomid) != TCL_OK) {
    return TCL_ERROR;
  }
  AtomIDList *aid = (AtomIDList *)clientData;
  aid->add(atomid-1);
  DebugM(4,"Atom ID " << atomid << " added to config list\n");
  return TCL_OK;
}


int ComputeGlobalMaster::Tcl_reconfig(ClientData clientData,
	Tcl_Interp *interp, int argc, char **) {
  if (argc != 1) {
    interp->result = "wrong # args";
    return TCL_ERROR;
  }
  int *reconfig = (int *)clientData;
  *reconfig = 1;
  DebugM(4,"Reconfiguration turned on\n");
  return TCL_OK;
}


int ComputeGlobalMaster::Tcl_loadcoords(ClientData clientData,
	Tcl_Interp *interp, int argc, char *argv[]) {
  if (argc != 2) {
    interp->result = "wrong # args";
    return TCL_ERROR;
  }
  char *vname = argv[1];
  ComputeGlobalMaster *self = (ComputeGlobalMaster *)clientData;
  char cmd[129];  int code;
  AtomIDList::iterator a_i = self->aid.begin();
  AtomIDList::iterator a_e = self->aid.end();
  PositionList::iterator p_i = self->p.begin();
  for ( ; a_i != a_e; ++a_i, ++p_i ) {
    sprintf(cmd, "set %s(%d) { %lg %lg %lg }", vname, (int)((*a_i)+1),
      (double)((*p_i).x),(double)((*p_i).y),(double)((*p_i).z));
    code = Tcl_Eval(interp,cmd);
    if (code != TCL_OK) {
      NAMD_die("TCL error in global force calculation!");
      return TCL_ERROR;
    }
  }
  return TCL_OK;
}


int ComputeGlobalMaster::Tcl_addforce(ClientData clientData,
	Tcl_Interp *interp, int argc, char *argv[]) {
  if (argc != 3) {
    interp->result = "wrong # args";
    return TCL_ERROR;
  }
  char **fstring;  int fnum;  int atomid;  double x, y, z;
  if (Tcl_GetInt(interp,argv[1],&atomid) != TCL_OK) {
    return TCL_ERROR;
  }
  if (Tcl_SplitList(interp, argv[2], &fnum, &fstring) != TCL_OK) {
    return TCL_ERROR;
  }
  if ( (fnum != 3) ||
       (Tcl_GetDouble(interp, fstring[0],&x) != TCL_OK) ||
       (Tcl_GetDouble(interp, fstring[1],&y) != TCL_OK) ||
       (Tcl_GetDouble(interp, fstring[2],&z) != TCL_OK) ) {
    interp->result = "force not a vector";
    free(fstring);
    return TCL_ERROR;
  }
  free(fstring);
  ComputeGlobalResultsMsg *msg = (ComputeGlobalResultsMsg *)clientData;
  msg->aid.add(atomid-1);
  msg->f.add(Vector(x,y,z));
  DebugM(4,"Atom ID " << atomid << " added to force list\n");
  return TCL_OK;
}


ComputeGlobalMaster::ComputeGlobalMaster(ComputeGlobal *h) {
  DebugM(3,"Constructing master\n");
  host = h;
  initialized = 0;
  msgcount = 0;
#ifdef NAMD_TCL
  interp = 0;
#endif
}

ComputeGlobalMaster::~ComputeGlobalMaster() {
  DebugM(3,"Destructing master\n");
#ifdef NAMD_TCL
  if ( interp ) Tcl_DeleteInterp(interp);
#endif
}

void ComputeGlobalMaster::recvData(ComputeGlobalDataMsg *msg) {
  DebugM(3,"Receiving data on master\n");
  // Check initialization and number of messages received
  if ( ! initialized )
  {
    delete msg;
    if ( ++msgcount == CNumPes() ) {
       msgcount = 0;
       initialize();
       initialized = 1;
    }
  }
  else {
    storedata(msg);
    if ( ++msgcount == CNumPes() ) {
       msgcount = 0;
       calculate();
       cleardata();
    }
  }
}

void ComputeGlobalMaster::initialize() {
  DebugM(4,"Initializing master\n");

  ComputeGlobalConfigMsg *msg =
	new (MsgIndex(ComputeGlobalConfigMsg)) ComputeGlobalConfigMsg;

#ifdef NAMD_TCL
  // Create interpreter
  interp = Tcl_CreateInterp();
  Tcl_CreateCommand(interp, "print", Tcl_print,
    (ClientData) NULL, (Tcl_CmdDeleteProc *) NULL);
  Tcl_CreateCommand(interp, "vecadd", proc_vecadd,
    (ClientData) NULL, (Tcl_CmdDeleteProc *) NULL);
  Tcl_CreateCommand(interp, "vecsub", proc_vecsub,
    (ClientData) NULL, (Tcl_CmdDeleteProc *) NULL);
  Tcl_CreateCommand(interp, "vecscale", proc_vecscale,
    (ClientData) NULL, (Tcl_CmdDeleteProc *) NULL);
/*
  Tcl_CreateCommand(interp, "transoffset", proc_transoffset,
    (ClientData) NULL, (Tcl_CmdDeleteProc *) NULL);
  Tcl_CreateCommand(interp, "transmult", proc_transmult,
    (ClientData) NULL, (Tcl_CmdDeleteProc *) NULL);
  Tcl_CreateCommand(interp, "vectrans", proc_vectrans,
    (ClientData) NULL, (Tcl_CmdDeleteProc *) NULL);
*/

  // Get the path for our script
  char *filename = Node::Object()->configList->find("globalForcesTcl")->data;

  // Call interpreter to determine requested atoms
  Tcl_CreateCommand(interp, "addatom", Tcl_addatom,
    (ClientData) &(msg->aid), (Tcl_CmdDeleteProc *) NULL);

  int code;
  code = Tcl_EvalFile(interp,filename);
  if (*interp->result != 0) CPrintf("TCL: %s\n",interp->result);
  if (code != TCL_OK) NAMD_die("TCL error in global force initialization!");

  Tcl_DeleteCommand(interp, "addatom");
#endif

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
  delete msg;
}

void ComputeGlobalMaster::cleardata() {
  aid.resize(0);
  p.resize(0);
}

void ComputeGlobalMaster::calculate() {
  DebugM(4,"Calculating forces on master\n");

  ComputeGlobalResultsMsg *msg =
	new (MsgIndex(ComputeGlobalResultsMsg)) ComputeGlobalResultsMsg;

#ifdef NAMD_TCL
  // Call interpreter to calculate forces
  Tcl_CreateCommand(interp, "loadcoords", Tcl_loadcoords,
    (ClientData) this, (Tcl_CmdDeleteProc *) NULL);
  Tcl_CreateCommand(interp, "addforce", Tcl_addforce,
    (ClientData) msg, (Tcl_CmdDeleteProc *) NULL);
  Tcl_CreateCommand(interp, "reconfig", Tcl_reconfig,
    (ClientData) &(msg->reconfig), (Tcl_CmdDeleteProc *) NULL);
  Tcl_CreateCommand(interp, "addatom", Tcl_addatom,
    (ClientData) &(msg->newaid), (Tcl_CmdDeleteProc *) NULL);

  char cmd[129];  int code;
  strcpy(cmd,"calcforces");  code = Tcl_Eval(interp,cmd);
  if (*interp->result != 0) CPrintf("TCL: %s\n",interp->result);
  if (code != TCL_OK) NAMD_die("TCL error in global force calculation!");

  Tcl_DeleteCommand(interp, "loadcoords");
  Tcl_DeleteCommand(interp, "addforce");
  Tcl_DeleteCommand(interp, "reconfig");
  Tcl_DeleteCommand(interp, "addatom");
#endif

  // Send results to clients
  DebugM(3,"Sending results (" << msg->aid.size() << " forces) on master\n");
  if ( msg->reconfig ) {
    DebugM(4,"Sending new configuration (" <<
			msg->newaid.size() << " atoms) on master\n");
  }
  host->comm->sendComputeGlobalResults(msg);
}


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
  master = ( CMyPe() ? 0 : new ComputeGlobalMaster(this) );
  comm = m;
  configured = 0;
}

ComputeGlobal::~ComputeGlobal()
{
  delete master;
}

void ComputeGlobal::recvConfig(ComputeGlobalConfigMsg *msg) {
  DebugM(4,"Receiving configuration (" <<
			msg->aid.size() << " atoms) on client\n");
  aid = msg->aid;
  delete msg;
  configured = 1;
  sendData();
}

void ComputeGlobal::recvResults(ComputeGlobalResultsMsg *msg) {
  DebugM(3,"Receiving results (" << msg->aid.size() << " forces) on client\n");

  // Get reconfiguration if present
  if ( msg->reconfig ) {
    DebugM(4,"Receiving new configuration (" <<
			msg->newaid.size() << " atoms) on client\n");
    aid = msg->newaid;
  }

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

  for (ap = ap.begin(); ap != ap.end(); ap++) {
    (*ap).forceBox->close(&((*ap).r));
  }
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
 *	$Revision: 1.1 $	$Date: 1997/12/19 23:48:45 $
 *
 ***************************************************************************
 * REVISION HISTORY:
 *
 * $Log: ComputeGlobal.C,v $
 * Revision 1.1  1997/12/19 23:48:45  jim
 * Added Tcl interface for calculating forces.
 *
 *
 ***************************************************************************/
