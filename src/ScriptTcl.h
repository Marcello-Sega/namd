/**
***  Copyright (c) 1995, 1996, 1997, 1998, 1999, 2000 by
***  The Board of Trustees of the University of Illinois.
***  All rights reserved.
**/

/*
   Modifies SimParameters settings during run.
*/

#ifndef SCRIPTTCL_H
#define SCRIPTTCL_H

#include "converse.h"
#include "NamdTypes.h"
#include "Broadcasts.h"

#ifdef NAMD_TCL
#include <tcl.h>
#endif

class ConfigList;

class ScriptTcl {
public:
  ScriptTcl();
  ~ScriptTcl();
  void awaken(void) { CthAwaken(thread); }
  void run(char *filename, ConfigList *configList);
private:
  char *scriptFile;
  ConfigList *config;
  CthThread thread;
  static void threadRun(ScriptTcl*);
  void suspend(void) { CthSuspend(); }
  void algorithm();
  int runWasCalled;
  void barrier();
  void runController(int task);
  void setParameter(const char* param, const char* value);
  void setParameter(const char* param, int value);
#ifdef NAMD_TCL
  friend class Controller;
  Tcl_Interp *interp;
  static int Tcl_print(ClientData, Tcl_Interp *, int, char **);
  static int Tcl_config(ClientData, Tcl_Interp *, int, char **);
  static int Tcl_param(ClientData, Tcl_Interp *, int, char **);
  static int Tcl_run(ClientData, Tcl_Interp *, int, char **);
  static int Tcl_move(ClientData, Tcl_Interp *, int, char **);
  static int Tcl_output(ClientData, Tcl_Interp *, int, char **);
  static int Tcl_callback(ClientData, Tcl_Interp *, int, char **);
  char *callbackname;
  void doCallback(const char *labels, const char *data);
  int doCallback() { return ! ! callbackname; }
#endif
};

#endif

