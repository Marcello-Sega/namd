/**
***  Copyright (c) 1995, 1996, 1997, 1998, 1999, 2000 by
***  The Board of Trustees of the University of Illinois.
***  All rights reserved.
**/

/*
   Forwards atoms to master node for force evaluation.
*/

#ifndef COMPUTETCL_H
#define COMPUTETCL_H

#include "ComputeGlobalMaster.h"
#include "NamdTypes.h"

class ComputeGlobalConfigMsg;
class ComputeGlobalDataMsg;
class ComputeGlobalResultsMsg;
class ComputeGlobalMaster;
class ComputeMgr;

#ifdef NAMD_TCL
#include <tcl.h>
#endif

class ComputeTcl : public ComputeGlobalMaster {
private:
  friend class ComputeGlobal;
  ComputeTcl(ComputeGlobal *);
  ~ComputeTcl();
  virtual void initialize();
  virtual void calculate();
#ifdef NAMD_TCL
  Tcl_Interp *interp;
  static int Tcl_print(ClientData, Tcl_Interp *, int, char **);
  static int Tcl_atomid(ClientData, Tcl_Interp *, int, char **);
  static int Tcl_addatom(ClientData, Tcl_Interp *, int, char **);
  static int Tcl_addgroup(ClientData, Tcl_Interp *, int, char **);
  static int Tcl_reconfig(ClientData, Tcl_Interp *, int, char **);
  static int Tcl_loadcoords(ClientData, Tcl_Interp *, int, char **);
  static int Tcl_loadmasses(ClientData, Tcl_Interp *, int, char **);
  static int Tcl_addforce(ClientData, Tcl_Interp *, int, char **);
#endif
};

#endif

