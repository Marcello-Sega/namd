/**
***  Copyright (c) 1995, 1996, 1997, 1998, 1999, 2000 by
***  The Board of Trustees of the University of Illinois.
***  All rights reserved.
**/

/*
   Modifies SimParameters settings during run.
*/

#include "BackEnd.h"
#include "ScriptTcl.h"
#include "Broadcasts.h"
#include "ConfigList.h"
#include "Node.h"
#include "NamdState.h"
#include "Controller.h"
#include "SimParameters.h"
#include "Thread.h"
#include "ProcessorPrivate.h"
#include "PatchMgr.h"
#include "Measure.h"
#include <stdio.h>
#include <ctype.h>  // for isspace
#ifndef WIN32
#include <strings.h>
#endif

#ifdef NAMD_TCL
#include <tcl.h>
#endif

//#define DEBUGM
#define MIN_DEBUG_LEVEL 4
#include "Debug.h"

void ScriptTcl::suspend() {
  BackEnd::suspend();
}

void ScriptTcl::barrier() {
  BackEnd::barrier();
}

void ScriptTcl::initcheck() {
  if ( runWasCalled == 0 ) {
#ifdef NAMD_TCL
    CkPrintf("TCL: Suspending until startup complete.\n");
    Tcl_CreateCommand(interp, "param", Tcl_param,
      (ClientData) this, (Tcl_CmdDeleteProc *) NULL);
    Tcl_CreateCommand(interp, "unknown", Tcl_param,
      (ClientData) this, (Tcl_CmdDeleteProc *) NULL);
#endif
    runWasCalled = 1;

    state->configListInit(config);
    Node::Object()->saveMolDataPointers(state);
    Node::messageStartUp();
    suspend();
  }
}

void ScriptTcl::runController(int task) {
  scriptBarrier.publish(barrierStep++,task);
  suspend();
}

void ScriptTcl::setParameter(const char* param, const char* value) {
  ScriptParamMsg *msg = new ScriptParamMsg;
  strncpy(msg->param,param,MAX_SCRIPT_PARAM_SIZE);
  strncpy(msg->value,value,MAX_SCRIPT_PARAM_SIZE);
  (CProxy_Node(CpvAccess(BOCclass_group).node)).scriptParam(msg);
  barrier();
}

void ScriptTcl::setParameter(const char* param, int value) {
  ScriptParamMsg *msg = new ScriptParamMsg;
  strncpy(msg->param,param,MAX_SCRIPT_PARAM_SIZE);
  sprintf(msg->value,"%d",value);
  (CProxy_Node(CpvAccess(BOCclass_group).node)).scriptParam(msg);
  barrier();
}

#ifdef NAMD_TCL

int ScriptTcl::Tcl_exit(ClientData clientData,
	Tcl_Interp *, int argc, char *argv[]) {
  ScriptTcl *script = (ScriptTcl *)clientData;
  script->runController(SCRIPT_END);
  BackEnd::exit();
  return TCL_OK;
}

int ScriptTcl::Tcl_abort(ClientData,
	Tcl_Interp *, int argc, char *argv[]) {
  Tcl_DString msg;
  Tcl_DStringInit(&msg);
  Tcl_DStringAppend(&msg,"TCL:",-1);
  for ( int i = 1; i < argc; ++i ) {
    Tcl_DStringAppend(&msg," ",-1);
    Tcl_DStringAppend(&msg,argv[i],-1);
  }
  NAMD_die(Tcl_DStringValue(&msg));
  Tcl_DStringFree(&msg);
  return TCL_OK;
}

int ScriptTcl::Tcl_print(ClientData,
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

int ScriptTcl::Tcl_config(ClientData clientData,
	Tcl_Interp *interp, int argc, char *argv[]) {

// Needs to handle the following cases as passed in by Tcl:
//    name data #comment
//    name=data #comment
//    name= data #comment
//    name =data #comment
//    name = data #comment
//    name data1 data2 data3 #comment
//    name=data1 data2 data3 #comment
//    name= data1 data2 data3 #comment
//    name =data1 data2 data3 #comment
//    name = data1 data2 data3 #comment
//    name { data1 data2 data3 } #comment
//    name { data1 data2 data3 } #comment
//    name { data1 data2 # data3 } #comment
//    name {data1 data2 # data3 } #comment
// Do not try to handle "data#comment" in any form.
// The '#' start of any comments will *always* be a new argv.
// The name will *always* be contained in argv[1].

  // allocate storage for data string
  int arglen = 1;  int ai;
  for (ai=1; ai<argc; ++ai) { arglen += strlen(argv[ai]) + 1; }
  char *data = new char[arglen];  *data = 0;

  // find the end of the name
  char *name, *s;
  name = argv[1];
  for ( s = name; *s && *s != '='; ++s );

  // eliminate any comment
  for (ai=2; ai<argc; ++ai) { if (argv[ai][0] == '#') argc = ai; }

  // concatenate all the data items
  ai = 2;
  if ( *s ) { *s = 0; ++s; strcat(data,s); }  // name=data or name=
  else if ( ai < argc && argv[ai][0] == '=' ) {  // name =data or name =
    strcat(data,argv[ai]+1);
    ++ai;
  }
  for ( ; ai<argc; ++ai) {
    if ( data[0] ) { strcat(data," "); }
    strcat(data,argv[ai]);
  }

  if ( ! *name || ! *data ) {
    delete [] data;
    Tcl_SetResult(interp,"error parsing config file",TCL_VOLATILE);
    return TCL_ERROR;
  }

  ScriptTcl *script = (ScriptTcl *)clientData;
  script->config->add_element( name, strlen(name), data, strlen(data) );

  delete [] data;
  return TCL_OK;
}

int ScriptTcl::Tcl_param(ClientData clientData,
        Tcl_Interp *interp, int argc, char *argv[]) {
  if (argc != 3 && argc != 5) {
    Tcl_SetResult(interp,"wrong # args",TCL_VOLATILE);
    return TCL_ERROR;
  }

  char *param = argv[1];
  if ( strlen(param) + 1 > MAX_SCRIPT_PARAM_SIZE ) {
    Tcl_SetResult(interp,"parameter name too long",TCL_VOLATILE);
    return TCL_ERROR;
  }

  char value[MAX_SCRIPT_PARAM_SIZE];
  int arglen = strlen(argv[2]) + 1;
  if ( argc == 5 ) arglen += strlen(argv[3]) + strlen(argv[4]) + 2;
  if ( arglen > MAX_SCRIPT_PARAM_SIZE ) {
    Tcl_SetResult(interp,"parameter value too long",TCL_VOLATILE);
    return TCL_ERROR;
  }
  if ( argc == 3 ) sprintf(value,"%s",argv[2]);
  if ( argc == 5 ) sprintf(value,"%s %s %s",argv[2],argv[3],argv[4]);

  iout << "TCL: Setting parameter " << param << " to " << value << "\n" << endi;

  ScriptTcl *script = (ScriptTcl *)clientData;
  script->setParameter(param,value);

  return TCL_OK;
}

int ScriptTcl::Tcl_reinitvels(ClientData clientData,
        Tcl_Interp *interp, int argc, char *argv[]) {
  ScriptTcl *script = (ScriptTcl *)clientData;
  script->initcheck();
  if (argc != 2) {
    Tcl_SetResult(interp,"wrong # args",TCL_VOLATILE);
    return TCL_ERROR;
  }
  char *temp = argv[1];

  script->setParameter("initialTemp",temp);

  script->runController(SCRIPT_REINITVELS);

  return TCL_OK;
}

int ScriptTcl::Tcl_run(ClientData clientData,
	Tcl_Interp *interp, int argc, char *argv[]) {
  ScriptTcl *script = (ScriptTcl *)clientData;
  script->initcheck();
  if (argc != 2) {
    Tcl_SetResult(interp,"wrong # args",TCL_VOLATILE);
    return TCL_ERROR;
  }
  int numsteps;
  if (Tcl_GetInt(interp,argv[1],&numsteps) != TCL_OK) {
    return TCL_ERROR;
  }
  if (numsteps < 0) {
    Tcl_SetResult(interp,"number of steps must be non-negative",TCL_VOLATILE);
    return TCL_ERROR;
  }
  SimParameters *simParams = Node::Object()->simParameters;
  if (numsteps % simParams->stepsPerCycle) {
    Tcl_SetResult(interp,"number of steps must be a multiple of stepsPerCycle",TCL_VOLATILE);
    return TCL_ERROR;
  }
  iout << "TCL: Running for " << numsteps << " steps\n" << endi;

  script->setParameter("numsteps",simParams->firstTimestep + numsteps);

  script->runController(SCRIPT_RUN);

  script->setParameter("firsttimestep",simParams->N);

  return TCL_OK;
}

int ScriptTcl::Tcl_minimize(ClientData clientData,
	Tcl_Interp *interp, int argc, char *argv[]) {
  ScriptTcl *script = (ScriptTcl *)clientData;
  script->initcheck();
  if (argc != 2) {
    Tcl_SetResult(interp,"wrong # args",TCL_VOLATILE);
    return TCL_ERROR;
  }
  int numsteps;
  if (Tcl_GetInt(interp,argv[1],&numsteps) != TCL_OK) {
    return TCL_ERROR;
  }
  if (numsteps < 0) {
    Tcl_SetResult(interp,"number of steps must be non-negative",TCL_VOLATILE);
    return TCL_ERROR;
  }
  SimParameters *simParams = Node::Object()->simParameters;
  if (numsteps % simParams->stepsPerCycle) {
    Tcl_SetResult(interp,"number of steps must be a multiple of stepsPerCycle",TCL_VOLATILE);
    return TCL_ERROR;
  }
  iout << "TCL: Minimizing for " << numsteps << " steps\n" << endi;

  script->setParameter("numsteps",simParams->firstTimestep + numsteps);

  script->runController(SCRIPT_MINIMIZE);

  script->setParameter("firsttimestep",simParams->N);

  return TCL_OK;
}

int ScriptTcl::Tcl_move(ClientData clientData,
	Tcl_Interp *interp, int argc, char *argv[]) {
  ScriptTcl *script = (ScriptTcl *)clientData;
  script->initcheck();
  if (argc != 4) {
    Tcl_SetResult(interp,"wrong # args",TCL_VOLATILE);
    return TCL_ERROR;
  }
  char **fstring;  int fnum;  int atomid;  int moveto;  double x, y, z;
  if (Tcl_GetInt(interp,argv[1],&atomid) != TCL_OK) return TCL_ERROR;
  if (argv[2][0]=='t' && argv[2][1]=='o' && argv[2][2]==0) moveto = 1;
  else if (argv[2][0]=='b' && argv[2][1]=='y' && argv[2][2]==0) moveto = 0;
  else {
    Tcl_SetResult(interp,"syntax is 'move <id> to|by {<x> <y> <z>}'",TCL_VOLATILE);
    return TCL_ERROR;
  }
  if (Tcl_SplitList(interp, argv[3], &fnum, &fstring) != TCL_OK) {
    return TCL_ERROR;
  }
  if ( (fnum != 3) ||
       (Tcl_GetDouble(interp, fstring[0],&x) != TCL_OK) ||
       (Tcl_GetDouble(interp, fstring[1],&y) != TCL_OK) ||
       (Tcl_GetDouble(interp, fstring[2],&z) != TCL_OK) ) {
    Tcl_SetResult(interp,"third argument not a vector",TCL_VOLATILE);
    Tcl_Free((char*)fstring);
    return TCL_ERROR;
  }
  Tcl_Free((char*)fstring);

  SimParameters *simParams = Node::Object()->simParameters;

  iout << "TCL: Moving atom " << atomid << " ";
  if ( moveto ) iout << "to"; else iout << "by";
  iout << " " << Vector(x,y,z) << ".\n" << endi;

  MoveAtomMsg *msg = new MoveAtomMsg;
  msg->atomid = atomid - 1;
  msg->moveto = moveto;
  msg->coord = Vector(x,y,z);
  (CProxy_PatchMgr(CpvAccess(BOCclass_group).patchMgr)).moveAtom(msg);

  script->barrier();

  return TCL_OK;
}

int ScriptTcl::Tcl_output(ClientData clientData,
        Tcl_Interp *interp, int argc, char *argv[]) {
  ScriptTcl *script = (ScriptTcl *)clientData;
  script->initcheck();
  if (argc != 2) {
    Tcl_SetResult(interp,"wrong # args",TCL_VOLATILE);
    return TCL_ERROR;
  }
  if (strlen(argv[1]) > MAX_SCRIPT_PARAM_SIZE) {
    Tcl_SetResult(interp,"file name too long",TCL_VOLATILE);
    return TCL_ERROR;
  }

  SimParameters *simParams = Node::Object()->simParameters;

  char oldname[MAX_SCRIPT_PARAM_SIZE+1];
  strncpy(oldname,simParams->outputFilename,MAX_SCRIPT_PARAM_SIZE);

  script->setParameter("outputname",argv[1]);

  iout << "TCL: Writing to files with basename " <<
		simParams->outputFilename << ".\n" << endi;

  script->runController(SCRIPT_OUTPUT);

  script->setParameter("outputname",oldname);

  return TCL_OK;
}

void ScriptTcl::measure(Vector *c) {
  Measure::createCommands(interp);
  Node::Object()->coords = c;
  Tcl_Eval(interp,measure_command);
  Node::Object()->coords = 0;
  Measure::deleteCommands(interp);
}

int ScriptTcl::Tcl_measure(ClientData clientData,
        Tcl_Interp *interp, int argc, char *argv[]) {
  ScriptTcl *script = (ScriptTcl *)clientData;
  script->initcheck();
  if (argc != 2) {
    Tcl_SetResult(interp,"wrong # args",TCL_VOLATILE);
    return TCL_ERROR;
  }
  script->measure_command = argv[1];

  script->runController(SCRIPT_MEASURE);

  return TCL_OK;
}

int ScriptTcl::Tcl_checkpoint(ClientData clientData,
        Tcl_Interp *interp, int argc, char *argv[]) {
  ScriptTcl *script = (ScriptTcl *)clientData;
  script->initcheck();
  if (argc != 1) {
    Tcl_SetResult(interp,"wrong # args",TCL_VOLATILE);
    return TCL_ERROR;
  }

  script->runController(SCRIPT_CHECKPOINT);

  return TCL_OK;
}

int ScriptTcl::Tcl_revert(ClientData clientData,
        Tcl_Interp *interp, int argc, char *argv[]) {
  ScriptTcl *script = (ScriptTcl *)clientData;
  script->initcheck();
  if (argc != 1) {
    Tcl_SetResult(interp,"wrong # args",TCL_VOLATILE);
    return TCL_ERROR;
  }

  script->runController(SCRIPT_REVERT);

  return TCL_OK;
}

int ScriptTcl::Tcl_callback(ClientData clientData,
	Tcl_Interp *interp, int argc, char *argv[]) {
  ScriptTcl *script = (ScriptTcl *)clientData;
  if (argc != 2) {
    Tcl_SetResult(interp,"wrong # args",TCL_VOLATILE);
    return TCL_ERROR;
  }

  delete [] script->callbackname;
  script->callbackname = new char[strlen(argv[1])+1];
  strcpy(script->callbackname,argv[1]);

  iout << "TCL: Reduction callback proc set to " <<
			script->callbackname << "\n" << endi;

  return TCL_OK;
}

void ScriptTcl::doCallback(const char *labels, const char *data) {
  if ( ! callbackname ) return;
  int len = strlen(callbackname) + strlen(labels) + strlen(data) + 7;
  char *cmd = new char[len];
  sprintf(cmd, "%s {%s} {%s}", callbackname, labels, data);
  int rval = Tcl_Eval(interp,cmd);
  delete [] cmd;
  if (rval != TCL_OK) {
    char *errorInfo = Tcl_GetVar(interp,"errorInfo",0);
    NAMD_die(errorInfo);
  }
}

#endif  // NAMD_TCL


ScriptTcl::ScriptTcl() : scriptBarrier(scriptBarrierTag) {
  DebugM(3,"Constructing ScriptTcl\n");
#ifdef NAMD_TCL
  interp = 0;
  callbackname = 0;
#endif
  state = new NamdState;
  barrierStep = 0;
}

ScriptTcl::~ScriptTcl() {
  DebugM(3,"Destructing ScriptTcl\n");
#ifdef NAMD_TCL
  if ( interp ) Tcl_DeleteInterp(interp);
  delete [] callbackname;
#endif
}

void ScriptTcl::run(char *filename, ConfigList *)
{
    scriptFile = filename;
    algorithm();
}

void ScriptTcl::algorithm() {
  DebugM(4,"Running ScriptTcl\n");

  runWasCalled = 0;

#ifdef NAMD_TCL
  config = new ConfigList;

  // Create interpreter
  interp = Tcl_CreateInterp();
  Tcl_CreateCommand(interp, "exit", Tcl_exit,
    (ClientData) this, (Tcl_CmdDeleteProc *) NULL);
  Tcl_CreateCommand(interp, "abort", Tcl_abort,
    (ClientData) NULL, (Tcl_CmdDeleteProc *) NULL);
  Tcl_CreateCommand(interp, "print", Tcl_print,
    (ClientData) NULL, (Tcl_CmdDeleteProc *) NULL);
  Tcl_CreateCommand(interp, "unknown", Tcl_config,
    (ClientData) this, (Tcl_CmdDeleteProc *) NULL);
  Tcl_CreateCommand(interp, "param", Tcl_config,
    (ClientData) this, (Tcl_CmdDeleteProc *) NULL);
  Tcl_CreateCommand(interp, "run", Tcl_run,
    (ClientData) this, (Tcl_CmdDeleteProc *) NULL);
  Tcl_CreateCommand(interp, "minimize", Tcl_minimize,
    (ClientData) this, (Tcl_CmdDeleteProc *) NULL);
  Tcl_CreateCommand(interp, "move", Tcl_move,
    (ClientData) this, (Tcl_CmdDeleteProc *) NULL);
  Tcl_CreateCommand(interp, "output", Tcl_output,
    (ClientData) this, (Tcl_CmdDeleteProc *) NULL);
  Tcl_CreateCommand(interp, "measure", Tcl_measure,
    (ClientData) this, (Tcl_CmdDeleteProc *) NULL);
  Tcl_CreateCommand(interp, "checkpoint", Tcl_checkpoint,
    (ClientData) this, (Tcl_CmdDeleteProc *) NULL);
  Tcl_CreateCommand(interp, "revert", Tcl_revert,
    (ClientData) this, (Tcl_CmdDeleteProc *) NULL);
  Tcl_CreateCommand(interp, "reinitvels", Tcl_reinitvels,
    (ClientData) this, (Tcl_CmdDeleteProc *) NULL);
  Tcl_CreateCommand(interp, "callback", Tcl_callback,
    (ClientData) this, (Tcl_CmdDeleteProc *) NULL);

  int code = Tcl_EvalFile(interp,scriptFile);
  char *result = Tcl_GetStringResult(interp);
  if (*result != 0) CkPrintf("TCL: %s\n",result);
  if (code != TCL_OK) {
    char *errorInfo = Tcl_GetVar(interp,"errorInfo",0);
    NAMD_die(errorInfo);
  }

#else
  if ( NULL == scriptFile || NULL == (config = new ConfigList(scriptFile)) ) {
    NAMD_die("Simulation config file is empty.");
  }
#endif

  if (runWasCalled == 0) {
    initcheck();
    SimParameters *simParams = Node::Object()->simParameters;
    if ( simParams->minimizeCGOn ) runController(SCRIPT_MINIMIZE);
    else runController(SCRIPT_RUN);
  }

  runController(SCRIPT_END);

}

