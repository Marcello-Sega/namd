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

void ScriptTcl::runController(int task) {
  scriptBarrier.publish(barrierStep++,task);
  suspend();
}

void ScriptTcl::setParameter(const char* param, const char* value) {
  ScriptParamMsg *msg = new ScriptParamMsg;
  strncpy(msg->param,param,MAX_SCRIPT_PARAM_SIZE);
  strncpy(msg->value,value,MAX_SCRIPT_PARAM_SIZE);
  CProxy_Node(CpvAccess(BOCclass_group).node).scriptParam(msg);
  barrier();
}

void ScriptTcl::setParameter(const char* param, int value) {
  ScriptParamMsg *msg = new ScriptParamMsg;
  strncpy(msg->param,param,MAX_SCRIPT_PARAM_SIZE);
  sprintf(msg->value,"%d",value);
  CProxy_Node(CpvAccess(BOCclass_group).node).scriptParam(msg);
  barrier();
}

#ifdef NAMD_TCL

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
  char *buf = Tcl_Merge(argc-1,argv+1);
  char *namestart, *nameend, *datastart, *dataend, *s;
  namestart = nameend = datastart = dataend = NULL;
  int spacecount = 0;
  int inbraces = 0;

    for (s = buf; *s; s++) {    // get to the end of the line
       if (*s == '#')                       // found a comment, so break
          { *s = 0; break; }
       if ( !isspace(*s) )    // dataend will always be the last non-blank char
          dataend = s;
       if ( !isspace(*s) && !namestart)     // found first character of name
          {namestart = s; continue; }
       if ( (isspace(*s)  || *s == '=') &&  // found last character of name
                 namestart && !nameend)
          nameend = s - 1;
       if ( !isspace(*s) && !datastart &&   // found the next char. after name
                 nameend) {
          if (*s == '=' && spacecount == 0) // an equals is allowed
             {spacecount++; continue; }     // but only once
          else if (*s == '{') {
            datastart = s;
            int escape_next = 0;
            int open_brace_count = 0;
            for(; *s; s++) {
              if (escape_next) { escape_next=0; continue; }
              if (*s == '\\' && ! escape_next) { escape_next=1; }
              if (*s == '{' && ! escape_next) { open_brace_count++; }
              if (*s == '}' && ! escape_next) { open_brace_count--; }
              if (! open_brace_count) { dataend = s; break; }
            }
            continue;
          }
          else
             {datastart = s; continue; }    // otherwise, use it
       }
    }

    if (!namestart || !nameend || !datastart || !dataend) {
      free(buf);
      interp->result = "error parsing config file";
      return TCL_ERROR;
    }

  ScriptTcl *script = (ScriptTcl *)clientData;
  script->config->add_element( namestart, nameend - namestart + 1,
                               datastart, dataend - datastart + 1 );

  free(buf);
  return TCL_OK;
}

int ScriptTcl::Tcl_param(ClientData clientData,
        Tcl_Interp *interp, int argc, char *argv[]) {
  if (argc != 3) {
    interp->result = "wrong # args";
    return TCL_ERROR;
  }
  char *param = argv[1];
  char *value = argv[2];

  iout << "TCL: Setting parameter " << param << " to " << value << "\n" << endi;

  ScriptTcl *script = (ScriptTcl *)clientData;
  script->setParameter(param,value);

  return TCL_OK;
}

int ScriptTcl::Tcl_run(ClientData clientData,
	Tcl_Interp *interp, int argc, char *argv[]) {
  ScriptTcl *script = (ScriptTcl *)clientData;
  if ( script->runWasCalled == 0 ) {
    CkPrintf("TCL: Run called, suspending until startup complete.\n");
    Tcl_CreateCommand(interp, "param", Tcl_param,
      (ClientData) script, (Tcl_CmdDeleteProc *) NULL);
    Tcl_CreateCommand(interp, "unknown", Tcl_param,
      (ClientData) script, (Tcl_CmdDeleteProc *) NULL);
    script->config->add_element("tcl",3,"on",2);
    script->runWasCalled = 1;

    script->state->configListInit(script->config);
    Node::Object()->saveMolDataPointers(script->state);
    Node::messageStartUp();
    script->suspend();
  }
  if (argc != 2) {
    interp->result = "wrong # args";
    return TCL_ERROR;
  }
  int numsteps;
  if (Tcl_GetInt(interp,argv[1],&numsteps) != TCL_OK) {
    return TCL_ERROR;
  }
  if (numsteps < 0) {
    interp->result = "number of steps must be non-negative";
    return TCL_ERROR;
  }
  SimParameters *simParams = Node::Object()->simParameters;
  if (numsteps % simParams->stepsPerCycle) {
    interp->result = "number of steps must be a multiple of stepsPerCycle";
    return TCL_ERROR;
  }
  iout << "TCL: Running for " << numsteps << " steps\n" << endi;

  script->setParameter("numsteps",simParams->firstTimestep + numsteps);

  script->runController(SCRIPT_RUN);

  script->setParameter("firsttimestep",simParams->N);

  return TCL_OK;
}

int ScriptTcl::Tcl_move(ClientData clientData,
	Tcl_Interp *interp, int argc, char *argv[]) {
  ScriptTcl *script = (ScriptTcl *)clientData;
  if (! script->runWasCalled) {
    interp->result = "called before run";
    return TCL_ERROR;
  }
  if (argc != 4) {
    interp->result = "wrong # args";
    return TCL_ERROR;
  }
  char **fstring;  int fnum;  int atomid;  int moveto;  double x, y, z;
  if (Tcl_GetInt(interp,argv[1],&atomid) != TCL_OK) return TCL_ERROR;
  if (argv[2][0]=='t' && argv[2][1]=='o' && argv[2][2]==0) moveto = 1;
  else if (argv[2][0]=='b' && argv[2][1]=='y' && argv[2][2]==0) moveto = 0;
  else {
    interp->result = "syntax is 'move <id> to|by {<x> <y> <z>}'";
    return TCL_ERROR;
  }
  if (Tcl_SplitList(interp, argv[3], &fnum, &fstring) != TCL_OK) {
    return TCL_ERROR;
  }
  if ( (fnum != 3) ||
       (Tcl_GetDouble(interp, fstring[0],&x) != TCL_OK) ||
       (Tcl_GetDouble(interp, fstring[1],&y) != TCL_OK) ||
       (Tcl_GetDouble(interp, fstring[2],&z) != TCL_OK) ) {
    interp->result = "third argument not a vector";
    free(fstring);
    return TCL_ERROR;
  }
  free(fstring);

  SimParameters *simParams = Node::Object()->simParameters;

  iout << "TCL: Moving atom " << atomid << " ";
  if ( moveto ) iout << "to"; else iout << "by";
  iout << " " << Vector(x,y,z) << ".\n" << endi;

  MoveAtomMsg *msg = new MoveAtomMsg;
  msg->atomid = atomid - 1;
  msg->moveto = moveto;
  msg->coord = Vector(x,y,z);
  CProxy_PatchMgr(CpvAccess(BOCclass_group).patchMgr).moveAtom(msg);

  script->barrier();

  return TCL_OK;
}

int ScriptTcl::Tcl_output(ClientData clientData,
        Tcl_Interp *interp, int argc, char *argv[]) {
  ScriptTcl *script = (ScriptTcl *)clientData;
  if (! script->runWasCalled) {
    interp->result = "called before run";
    return TCL_ERROR;
  }
  if (argc != 2) {
    interp->result = "wrong # args";
    return TCL_ERROR;
  }
  if (strlen(argv[1]) > MAX_SCRIPT_PARAM_SIZE) {
    interp->result = "file name too long";
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
  if (! script->runWasCalled) {
    interp->result = "called before run";
    return TCL_ERROR;
  }
  if (argc != 2) {
    interp->result = "wrong # args";
    return TCL_ERROR;
  }
  script->measure_command = argv[1];

  script->runController(SCRIPT_MEASURE);

  return TCL_OK;
}

int ScriptTcl::Tcl_callback(ClientData clientData,
	Tcl_Interp *interp, int argc, char *argv[]) {
  ScriptTcl *script = (ScriptTcl *)clientData;
  if (argc != 2) {
    interp->result = "wrong # args";
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
     char *errmsg = new char[strlen(interp->result) + 20];
     sprintf(errmsg,"Tcl callback: %s",interp->result);
     NAMD_die(errmsg);
     delete [] errmsg;
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

void ScriptTcl::run(char *filename, ConfigList *configList)
{
    scriptFile = filename;
    config = configList;

    algorithm();
}

void ScriptTcl::algorithm() {
  DebugM(4,"Running ScriptTcl\n");

#ifdef NAMD_TCL
  // Create interpreter
  interp = Tcl_CreateInterp();
//  if (Tcl_Init(interp) == TCL_ERROR) {
//    CkPrintf("Tcl startup error: %s\n", interp->result);
//  }
  Tcl_CreateCommand(interp, "print", Tcl_print,
    (ClientData) NULL, (Tcl_CmdDeleteProc *) NULL);
  Tcl_CreateCommand(interp, "unknown", Tcl_config,
    (ClientData) this, (Tcl_CmdDeleteProc *) NULL);
  Tcl_CreateCommand(interp, "param", Tcl_config,
    (ClientData) this, (Tcl_CmdDeleteProc *) NULL);
  Tcl_CreateCommand(interp, "run", Tcl_run,
    (ClientData) this, (Tcl_CmdDeleteProc *) NULL);
  Tcl_CreateCommand(interp, "move", Tcl_move,
    (ClientData) this, (Tcl_CmdDeleteProc *) NULL);
  Tcl_CreateCommand(interp, "output", Tcl_output,
    (ClientData) this, (Tcl_CmdDeleteProc *) NULL);
  Tcl_CreateCommand(interp, "measure", Tcl_measure,
    (ClientData) this, (Tcl_CmdDeleteProc *) NULL);
  Tcl_CreateCommand(interp, "callback", Tcl_callback,
    (ClientData) this, (Tcl_CmdDeleteProc *) NULL);

/*
  // Get the script
  StringList *script = Node::Object()->configList->find("tclScript");

  for ( ; script; script = script->next ) {
    int code;
    if ( script->data[0] == '{' ) code = Tcl_Eval(interp,script->data+1);
    else code = Tcl_EvalFile(interp,script->data);
    if (*interp->result != 0) CkPrintf("TCL: %s\n",interp->result);
    if (code != TCL_OK) NAMD_die("TCL error in script!");
  }
*/

  runWasCalled = 0;
  int code = Tcl_EvalFile(interp,scriptFile);
  if (*interp->result != 0) CkPrintf("TCL: %s\n",interp->result);
  if (code != TCL_OK) NAMD_die("TCL error in script!");
  if (runWasCalled == 0) {
    if (callbackname == 0) {
      CkPrintf("TCL: Exiting after processing config file.\n");
      Tcl_DeleteInterp(interp);
      interp = 0;
    }
    state->configListInit(config);
    Node::Object()->saveMolDataPointers(state);
    Node::messageStartUp();
    suspend();
    return;
  } else {
    runController(SCRIPT_END);
  }

#else

  NAMD_die("Sorry, TCL scripting not available; built without TCL.");

#endif

}

