/**
***  Copyright (c) 1995, 1996, 1997, 1998, 1999, 2000 by
***  The Board of Trustees of the University of Illinois.
***  All rights reserved.
**/

#include "ComputeTclBC.h"
#include "Node.h"
#include "SimParameters.h"
#include "Patch.h"

#ifdef NAMD_TCL
#include <tcl.h>
#endif
#include "TclCommands.h"


ComputeTclBC::ComputeTclBC(ComputeID c)
  : ComputeHomePatches(c), ap(patchList) {

  reduction = ReductionMgr::Object()->willSubmit(REDUCTIONS_BASIC);
  SimParameters *simParams = Node::Object()->simParameters;

#ifdef NAMD_TCL
  interp = Tcl_CreateInterp();
  Tcl_CreateCommand(interp, "print", Tcl_print,
    (ClientData) NULL, (Tcl_CmdDeleteProc *) NULL);
  Tcl_CreateCommand(interp, "vecadd", proc_vecadd,
    (ClientData) NULL, (Tcl_CmdDeleteProc *) NULL);
  Tcl_CreateCommand(interp, "vecsub", proc_vecsub,
    (ClientData) NULL, (Tcl_CmdDeleteProc *) NULL);
  Tcl_CreateCommand(interp, "vecscale", proc_vecscale,
    (ClientData) NULL, (Tcl_CmdDeleteProc *) NULL);
  Tcl_CreateCommand(interp, "getbond", proc_getbond,
    (ClientData) NULL, (Tcl_CmdDeleteProc *) NULL);
  Tcl_CreateCommand(interp, "getangle", proc_getangle,
    (ClientData) NULL, (Tcl_CmdDeleteProc *) NULL);
  Tcl_CreateCommand(interp, "getdihedral", proc_getdihedral,
    (ClientData) NULL, (Tcl_CmdDeleteProc *) NULL);
  Tcl_CreateCommand(interp, "anglegrad", proc_anglegrad,
    (ClientData) NULL, (Tcl_CmdDeleteProc *) NULL);
  Tcl_CreateCommand(interp, "dihedralgrad", proc_dihedralgrad,
    (ClientData) NULL, (Tcl_CmdDeleteProc *) NULL);

  // run script to define calcforces, etc.
  if ( simParams->tclBCScript ) {
    int code = Tcl_Eval(interp,simParams->tclBCScript);
    char *result = Tcl_GetStringResult(interp);
    if (result && *result != 0) CkPrintf("TCL: %s\n",result);
    if (code != TCL_OK) {
      char *errorInfo = Tcl_GetVar(interp,"errorInfo",0);
      NAMD_die(errorInfo);
    }
  } else NAMD_bug("tclBCScript pointer was NULL");

  // don't want these available until calcforces call
  Tcl_CreateObjCommand(interp, "nextatom", Tcl_nextatom,
    (ClientData) this, (Tcl_CmdDeleteProc *) NULL);
  Tcl_CreateObjCommand(interp, "getcoord", Tcl_getcoord,
    (ClientData) this, (Tcl_CmdDeleteProc *) NULL);
  Tcl_CreateObjCommand(interp, "getmass", Tcl_getmass,
    (ClientData) this, (Tcl_CmdDeleteProc *) NULL);
  Tcl_CreateObjCommand(interp, "getcharge", Tcl_getcharge,
    (ClientData) this, (Tcl_CmdDeleteProc *) NULL);
  Tcl_CreateObjCommand(interp, "getid", Tcl_getid,
    (ClientData) this, (Tcl_CmdDeleteProc *) NULL);
  Tcl_CreateObjCommand(interp, "addforce", Tcl_addforce,
    (ClientData) this, (Tcl_CmdDeleteProc *) NULL);
  Tcl_CreateObjCommand(interp, "addenergy", Tcl_addenergy,
    (ClientData) this, (Tcl_CmdDeleteProc *) NULL);

#else

  NAMD_die("Sorry, tclBC is not available; built without TCL.");

#endif

}


ComputeTclBC::~ComputeTclBC() {
#ifdef NAMD_TCL
  Tcl_DeleteInterp(interp);
#endif
  delete reduction;
}


void ComputeTclBC::doWork() {

  SimParameters *simParams = Node::Object()->simParameters;
  const Lattice & lattice = patchList[0].p->lattice;
  const int step = patchList[0].p->flags.step;
  char cmd[128];

  energy = 0;
  n_atom = -1;  // set initial flags for iteration by nextatom

#ifdef NAMD_TCL
  sprintf(cmd,"calcforces %d %d %s",step,hasPatchZero,simParams->tclBCArgs);
  int code = Tcl_Eval(interp,cmd);
  if (code != TCL_OK) {
    char *errorInfo = Tcl_GetVar(interp,"errorInfo",0);
    NAMD_die(errorInfo);
  }
  if (n_atom != -2) {
    NAMD_die("tclBCScript failed to call nextatom until failure");
  }
#endif

  reduction->item(REDUCTION_BC_ENERGY) += energy;
  reduction->submit();

}

#ifdef NAMD_TCL

int ComputeTclBC::Tcl_print(ClientData,
        Tcl_Interp *, int argc, char *argv[]) {
  Tcl_DString msg;
  Tcl_DStringInit(&msg);
  for ( int i = 1; i < argc; ++i ) {
    Tcl_DStringAppend(&msg," ",-1);
    Tcl_DStringAppend(&msg,argv[i],-1);
  }
  CkPrintf("TCL:%s\n",Tcl_DStringValue(&msg));
  Tcl_DStringFree(&msg);
  return TCL_OK;
}

int ComputeTclBC::Tcl_nextatom(ClientData clientData,
        Tcl_Interp *interp, int objc, Tcl_Obj * const objv[]) {
  if (objc != 1) {
    Tcl_SetResult(interp,"wrong # args",TCL_VOLATILE);
    return TCL_ERROR;
  }
  ComputeTclBC *self = (ComputeTclBC *)clientData;

  // n_atom = -2 after all atoms processed
  if (self->n_atom < -1) {
    Tcl_SetObjResult(interp, Tcl_NewIntObj((long)(0)));
    return TCL_OK;
  }

  // assume n_atom = -1 before first call
  while ( self->n_atom < 0 || ++self->i_atom >= self->n_atom ) {
    if ( self->n_atom < 0 ) {  // first call
      self->ap = self->ap.begin();
    } else {
      (*(self->ap)).positionBox->close(&(self->atoms));
      (*(self->ap)).forceBox->close(&((*(self->ap)).r));
      self->ap++;
    }
    if ( self->ap == self->ap.end() ) {
      self->n_atom = -2;  // set error condition
      Tcl_SetObjResult(interp, Tcl_NewIntObj((long)(0)));
      return TCL_OK;
    }
    self->i_atom = -1;
    self->n_atom = (*(self->ap)).p->getNumAtoms();
    self->fullatoms = (*(self->ap)).p->getAtomList().begin();
    self->atoms = (*(self->ap)).positionBox->open();
    (*(self->ap)).r = (*(self->ap)).forceBox->open();
    self->forces = (*(self->ap)).r->f[Results::normal];
  }

  Tcl_SetObjResult(interp, Tcl_NewIntObj((long)(1)));
  return TCL_OK;
}

int ComputeTclBC::Tcl_getcoord(ClientData clientData,
        Tcl_Interp *interp, int objc, Tcl_Obj * const objv[]) {
  if (objc != 1) {
    Tcl_SetResult(interp,"wrong # args",TCL_VOLATILE);
    return TCL_ERROR;
  }
  ComputeTclBC *self = (ComputeTclBC *)clientData;
  if ( self->n_atom <= 0 ) {
    Tcl_SetResult(interp,"no atom available",TCL_VOLATILE);
    return TCL_ERROR;
  }
  Tcl_Obj *newlist = Tcl_NewListObj(0, NULL);

  int i = self->i_atom;
  Tcl_ListObjAppendElement(interp, newlist,
      Tcl_NewDoubleObj((double)(self->atoms[i].position.x)));
  Tcl_ListObjAppendElement(interp, newlist,
      Tcl_NewDoubleObj((double)(self->atoms[i].position.y)));
  Tcl_ListObjAppendElement(interp, newlist,
      Tcl_NewDoubleObj((double)(self->atoms[i].position.z)));

  Tcl_SetObjResult(interp, newlist);
  return TCL_OK;
}

int ComputeTclBC::Tcl_getmass(ClientData clientData,
        Tcl_Interp *interp, int objc, Tcl_Obj * const objv[]) {
  if (objc != 1) {
    Tcl_SetResult(interp,"wrong # args",TCL_VOLATILE);
    return TCL_ERROR;
  }
  ComputeTclBC *self = (ComputeTclBC *)clientData;
  if ( self->n_atom <= 0 ) {
    Tcl_SetResult(interp,"no atom available",TCL_VOLATILE);
    return TCL_ERROR;
  }

  int i = self->i_atom;
  Tcl_SetObjResult(interp, Tcl_NewDoubleObj((double)(self->fullatoms[i].mass)));
  return TCL_OK;
}

int ComputeTclBC::Tcl_getcharge(ClientData clientData,
        Tcl_Interp *interp, int objc, Tcl_Obj * const objv[]) {
  if (objc != 1) {
    Tcl_SetResult(interp,"wrong # args",TCL_VOLATILE);
    return TCL_ERROR;
  }
  ComputeTclBC *self = (ComputeTclBC *)clientData;
  if ( self->n_atom <= 0 ) {
    Tcl_SetResult(interp,"no atom available",TCL_VOLATILE);
    return TCL_ERROR;
  }

  int i = self->i_atom;
  Tcl_SetObjResult(interp, Tcl_NewDoubleObj((double)(self->atoms[i].charge)));
  return TCL_OK;
}

int ComputeTclBC::Tcl_getid(ClientData clientData,
        Tcl_Interp *interp, int objc, Tcl_Obj * const objv[]) {
  if (objc != 1) {
    Tcl_SetResult(interp,"wrong # args",TCL_VOLATILE);
    return TCL_ERROR;
  }
  ComputeTclBC *self = (ComputeTclBC *)clientData;
  if ( self->n_atom <= 0 ) {
    Tcl_SetResult(interp,"no atom available",TCL_VOLATILE);
    return TCL_ERROR;
  }

  int i = self->i_atom;
  Tcl_SetObjResult(interp, Tcl_NewIntObj((long)(self->atoms[i].id)));
  return TCL_OK;
}

int ComputeTclBC::Tcl_addforce(ClientData clientData,
        Tcl_Interp *interp, int objc, Tcl_Obj * const objv[]) {
  if (objc != 2) {
    Tcl_SetResult(interp,"wrong # args",TCL_VOLATILE);
    return TCL_ERROR;
  }

  Tcl_Obj **force;  int fnum;  double x,y,z;
  if (Tcl_ListObjGetElements(interp, objv[1], &fnum, &force) != TCL_OK) {
    return TCL_ERROR;
  }
  if ( (fnum != 3) ||
       (Tcl_GetDoubleFromObj(interp, force[0],&x) != TCL_OK) ||
       (Tcl_GetDoubleFromObj(interp, force[1],&y) != TCL_OK) ||
       (Tcl_GetDoubleFromObj(interp, force[2],&z) != TCL_OK) ) {
    Tcl_SetResult(interp,"force not a vector",TCL_VOLATILE);
    return TCL_ERROR;
  }

  ComputeTclBC *self = (ComputeTclBC *)clientData;
  if ( self->n_atom <= 0 ) {
    Tcl_SetResult(interp,"no atom available",TCL_VOLATILE);
    return TCL_ERROR;
  }
  int i = self->i_atom;
  self->forces[i].x += x;
  self->forces[i].y += y;
  self->forces[i].z += z;

  return TCL_OK;
}

int ComputeTclBC::Tcl_addenergy(ClientData clientData,
        Tcl_Interp *interp, int objc, Tcl_Obj * const objv[]) {
  if (objc != 2) {
    Tcl_SetResult(interp,"wrong # args",TCL_VOLATILE);
    return TCL_ERROR;
  }

  double energy;
  if ( Tcl_GetDoubleFromObj(interp, objv[1], &energy) != TCL_OK ) {
    Tcl_SetResult(interp,"energy not a number",TCL_VOLATILE);
    return TCL_ERROR;
  }

  ComputeTclBC *self = (ComputeTclBC *)clientData;
  self->energy += energy;

  return TCL_OK;
}

#endif

