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
class NamdState;

class ScriptTcl {
public:
  ScriptTcl();
  ~ScriptTcl();
  void run(char *filename, ConfigList *configList);
  void measure(Vector *);
private:
  char *scriptFile;
  ConfigList *config;
  NamdState *state;
  void suspend(void);
  void algorithm();
  int runWasCalled;
  void barrier();
  void initcheck();
  SimpleBroadcastObject<int> scriptBarrier;
  int barrierStep;
  void runController(int task);
  void setParameter(const char* param, const char* value);
  void setParameter(const char* param, int value);
#ifdef NAMD_TCL
  friend class Controller;
  friend class GlobalMasterTcl;
  Tcl_Interp *interp;
  static int Tcl_exit(ClientData, Tcl_Interp *, int, char **);
  static int Tcl_abort(ClientData, Tcl_Interp *, int, char **);
  static int Tcl_print(ClientData, Tcl_Interp *, int, char **);
  static int Tcl_config(ClientData, Tcl_Interp *, int, char **);
  static int Tcl_param(ClientData, Tcl_Interp *, int, char **);
  static int Tcl_reinitvels(ClientData, Tcl_Interp *, int, char **);
  static int Tcl_run(ClientData, Tcl_Interp *, int, char **);
  static int Tcl_minimize(ClientData, Tcl_Interp *, int, char **);
  static int Tcl_move(ClientData, Tcl_Interp *, int, char **);
  static int Tcl_output(ClientData, Tcl_Interp *, int, char **);
  static int Tcl_measure(ClientData, Tcl_Interp *, int, char **);
  static int Tcl_checkpoint(ClientData, Tcl_Interp *, int, char **);
  static int Tcl_revert(ClientData, Tcl_Interp *, int, char **);
  static int Tcl_callback(ClientData, Tcl_Interp *, int, char **);
  char *callbackname;
  void doCallback(const char *labels, const char *data);
  int doCallback() { return ! ! callbackname; }
  char *measure_command;
#endif
};

#endif

