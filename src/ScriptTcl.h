/***************************************************************************/
/*      (C) Copyright 1996,1997 The Board of Trustees of the               */
/*                          University of Illinois                         */
/*                           All Rights Reserved                           */
/***************************************************************************/
/***************************************************************************
 * DESCRIPTION:  Modifies SimParameters settings during run.
 *
 ***************************************************************************/

#ifndef SCRIPTTCL_H
#define SCRIPTTCL_H

#include "converse.h"
#include "NamdTypes.h"
#include "Broadcasts.h"

#ifdef NAMD_TCL
#include <tcl.h>
#endif

class ScriptTcl {
public:
  ScriptTcl();
  ~ScriptTcl();
  void awaken(void) { CthAwaken(thread); }
  void run();
private:
  CthThread thread;
  static void threadRun(ScriptTcl*);
  void suspend(void) { CthSuspend(); }
  void algorithm();
  SimpleBroadcastObject<int> scriptBarrier;
  int barrierStep;
#ifdef NAMD_TCL
  Tcl_Interp *interp;
  static int Tcl_print(ClientData, Tcl_Interp *, int, char **);
  static int Tcl_param(ClientData, Tcl_Interp *, int, char **);
  static int Tcl_run(ClientData, Tcl_Interp *, int, char **);
  static int Tcl_output(ClientData, Tcl_Interp *, int, char **);
#endif
};

#endif
/***************************************************************************
 * RCS INFORMATION:
 *
 *	$RCSfile: ScriptTcl.h,v $
 *	$Author: jim $	$Locker:  $		$State: Exp $
 *	$Revision: 1.2 $	$Date: 1999/06/21 16:15:36 $
 *
 ***************************************************************************
 * REVISION HISTORY:
 *
 * $Log: ScriptTcl.h,v $
 * Revision 1.2  1999/06/21 16:15:36  jim
 * Improved scripting, run now ends and generates output.
 *
 * Revision 1.1  1999/05/27 18:38:58  jim
 * Files to implement general Tcl scripting.
 *
 *
 *
 ***************************************************************************/

