/***************************************************************************/
/*      (C) Copyright 1996,1997 The Board of Trustees of the               */
/*                          University of Illinois                         */
/*                           All Rights Reserved                           */
/***************************************************************************/
/***************************************************************************
 * DESCRIPTION:  Forwards atoms to master node for force evaluation.
 *
 ***************************************************************************/

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
  static int Tcl_addatom(ClientData, Tcl_Interp *, int, char **);
  static int Tcl_reconfig(ClientData, Tcl_Interp *, int, char **);
  static int Tcl_loadcoords(ClientData, Tcl_Interp *, int, char **);
  static int Tcl_loadmasses(ClientData, Tcl_Interp *, int, char **);
  static int Tcl_addforce(ClientData, Tcl_Interp *, int, char **);
#endif
};

#endif
/***************************************************************************
 * RCS INFORMATION:
 *
 *	$RCSfile: ComputeTcl.h,v $
 *	$Author: jim $	$Locker:  $		$State: Exp $
 *	$Revision: 1.1 $	$Date: 1998/02/10 05:35:05 $
 *
 ***************************************************************************
 * REVISION HISTORY:
 *
 * $Log: ComputeTcl.h,v $
 * Revision 1.1  1998/02/10 05:35:05  jim
 * Split ComputeGlobal into different classes and files.
 * Switched globalForces and globalForcesTcl to tclForces and tclForcesScript.
 * Added (soon to be used) freeEnergy and freeEnergyConfig.
 *
 *
 *
 ***************************************************************************/

