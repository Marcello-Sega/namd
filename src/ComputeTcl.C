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
#include "ComputeTcl.h"
#include "ComputeGlobal.h"
#include "ComputeGlobalMsgs.h"
#include "Molecule.h"
#include "ReductionMgr.h"
#include "ComputeMgr.h"
#include "ComputeMgr.top.h"
#include <stdio.h>

#ifdef NAMD_TCL
#include <tcl.h>
#include <tclExtend.h>
#include "TclCommands.h"
#endif

// #define DEBUGM
#define MIN_DEBUG_LEVEL 4
#include "Debug.h"


#ifdef NAMD_TCL
int ComputeTcl::Tcl_print(ClientData,
	Tcl_Interp *, int argc, char *argv[]) {
  char *msg = Tcl_Merge(argc-1,argv+1);
  CPrintf("TCL: %s\n",msg);
  free(msg);
  return TCL_OK;
}


int ComputeTcl::Tcl_atomid(ClientData clientData,
	Tcl_Interp *interp, int argc, char *argv[]) {
  if (argc != 4) {
    interp->result = "wrong # args";
    return TCL_ERROR;
  }
  char *segid = argv[1];
  int resid;
  if (Tcl_GetInt(interp,argv[2],&resid) != TCL_OK) {
    return TCL_ERROR;
  }
  char *aname = argv[3];

  Molecule *mol = (Molecule *)clientData;
  int atomid = mol->get_atom_from_name(segid,resid,aname);

  if (atomid < 0) {
    interp->result = "atom not found";
    return TCL_ERROR;
  }
  atomid += 1;

  char s[10];  sprintf(s,"%d",atomid);
  Tcl_SetResult(interp,s,TCL_VOLATILE);
  DebugM(4,"Atom ID " << atomid << " identified by name\n");
  return TCL_OK;
}


int ComputeTcl::Tcl_addatom(ClientData clientData,
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


int ComputeTcl::Tcl_reconfig(ClientData clientData,
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


int ComputeTcl::Tcl_loadcoords(ClientData clientData,
	Tcl_Interp *interp, int argc, char *argv[]) {
  if (argc != 2) {
    interp->result = "wrong # args";
    return TCL_ERROR;
  }
  char *vname = argv[1];
  ComputeTcl *self = (ComputeTcl *)clientData;
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


int ComputeTcl::Tcl_loadmasses(ClientData clientData,
	Tcl_Interp *interp, int argc, char *argv[]) {
  if (argc != 2) {
    interp->result = "wrong # args";
    return TCL_ERROR;
  }
  char *vname = argv[1];
  ComputeTcl *self = (ComputeTcl *)clientData;
  char cmd[129];  int code;
  Molecule *mol = Node::Object()->molecule;
  AtomIDList::iterator a_i = self->aid.begin();
  AtomIDList::iterator a_e = self->aid.end();
  for ( ; a_i != a_e; ++a_i) {
    sprintf(cmd, "set %s(%d) %lg", vname, (int)((*a_i)+1),
      (double)(mol->atommass(*a_i)) );
    code = Tcl_Eval(interp,cmd);
    if (code != TCL_OK) {
      NAMD_die("TCL error in global force calculation!");
      return TCL_ERROR;
    }
  }
  return TCL_OK;
}


int ComputeTcl::Tcl_addforce(ClientData clientData,
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
#endif


ComputeTcl::ComputeTcl(ComputeGlobal *h) : ComputeGlobalMaster(h) {
  DebugM(3,"Constructing ComputeTcl\n");
#ifdef NAMD_TCL
  interp = 0;
#endif
}

ComputeTcl::~ComputeTcl() {
  DebugM(3,"Destructing ComputeTcl\n");
#ifdef NAMD_TCL
  if ( interp ) Tcl_DeleteInterp(interp);
#endif
}


void ComputeTcl::initialize() {
  DebugM(4,"Initializing master\n");

  ComputeGlobalConfigMsg *msg =
	new (MsgIndex(ComputeGlobalConfigMsg)) ComputeGlobalConfigMsg;

#ifdef NAMD_TCL
  // Create interpreter
  interp = Tcl_CreateInterp();
  if (Tcl_Init(interp) == TCL_ERROR) {
    CPrintf("Tcl startup error: %\n", interp->result);
    }
  if (Tclx_Init(interp) == TCL_ERROR) {
    CPrintf("Tcl-X startup error: %s\n", interp->result);
    } else {
      Tcl_StaticPackage(interp, "Tclx", Tclx_Init, Tclx_SafeInit);
    }
  Tcl_CreateCommand(interp, "print", Tcl_print,
    (ClientData) NULL, (Tcl_CmdDeleteProc *) NULL);
  Tcl_CreateCommand(interp, "atomid", Tcl_atomid,
    (ClientData) (Node::Object()->molecule), (Tcl_CmdDeleteProc *) NULL);
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
  char *filename = Node::Object()->configList->find("tclForcesScript")->data;

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


void ComputeTcl::calculate() {
  DebugM(4,"Calculating forces on master\n");

  ComputeGlobalResultsMsg *msg =
	new (MsgIndex(ComputeGlobalResultsMsg)) ComputeGlobalResultsMsg;

#ifdef NAMD_TCL
  // Call interpreter to calculate forces
  Tcl_CreateCommand(interp, "loadcoords", Tcl_loadcoords,
    (ClientData) this, (Tcl_CmdDeleteProc *) NULL);
  Tcl_CreateCommand(interp, "loadmasses", Tcl_loadmasses,
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
  Tcl_DeleteCommand(interp, "loadmasses");
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


/***************************************************************************
 * RCS INFORMATION:
 *
 *	$RCSfile $
 *	$Author $	$Locker:  $		$State: Exp $
 *	$Revision: 1.3 $	$Date: 1998/02/11 09:13:25 $
 *
 ***************************************************************************
 * REVISION HISTORY:
 *
 * $Log: ComputeTcl.C,v $
 * Revision 1.3  1998/02/11 09:13:25  jim
 * Added atomid command to tclForces.  Finds id from segname, resid, atomname.
 *
 * Revision 1.2  1998/02/10 06:45:10  jim
 * Added class ComputeFreeEnergy.
 *
 * Revision 1.1  1998/02/10 05:35:04  jim
 * Split ComputeGlobal into different classes and files.
 * Switched globalForces and globalForcesTcl to tclForces and tclForcesScript.
 * Added (soon to be used) freeEnergy and freeEnergyConfig.
 *
 *
 *
 ***************************************************************************/
