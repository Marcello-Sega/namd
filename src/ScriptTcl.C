/***************************************************************************/
/*       (C) Copyright 1996,1997 The Board of Trustees of the              */
/*                          University of Illinois                         */
/*                           All Rights Reserved                           */
/***************************************************************************/
/***************************************************************************
 * DESCRIPTION:  Modifies SimParameters settings during run.
 *
 ***************************************************************************/

#include "ScriptTcl.h"
#include "Broadcasts.h"
#include "ConfigList.h"
#include "Node.h"
#include "SimParameters.h"
#include "Thread.h"
#include "ProcessorPrivate.h"
#include "PatchMgr.h"
#include <stdio.h>

#ifdef NAMD_TCL
#include <tcl.h>
#endif

//#define DEBUGM
#define MIN_DEBUG_LEVEL 4
#include "Debug.h"

#ifdef NAMD_TCL

int ScriptTcl::Tcl_print(ClientData,
	Tcl_Interp *, int argc, char *argv[]) {
  char *msg = Tcl_Merge(argc-1,argv+1);
  CkPrintf("TCL: %s\n",msg);
  free(msg);
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

  ScriptParamMsg *msg = new ScriptParamMsg;
  strncpy(msg->param,param,MAX_SCRIPT_PARAM_SIZE);
  strncpy(msg->value,value,MAX_SCRIPT_PARAM_SIZE);
  CProxy_Node(CpvAccess(BOCclass_group).node).scriptParam(msg);

  Node::Object()->enableScriptBarrier();
  ScriptTcl *script = (ScriptTcl *)clientData;
  script->suspend();
  return TCL_OK;
}

int ScriptTcl::Tcl_run(ClientData clientData,
	Tcl_Interp *interp, int argc, char *argv[]) {
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
  ScriptTcl *script = (ScriptTcl *)clientData;
  if (numsteps % simParams->stepsPerCycle) {
    interp->result = "number of steps must be a multiple of stepsPerCycle";
    return TCL_ERROR;
  }
  iout << "TCL: Running for " << numsteps << " steps\n" << endi;

  ScriptParamMsg *msg = new ScriptParamMsg;
  sprintf(msg->param,"numsteps");
  sprintf(msg->value,"%d",simParams->firstTimestep + numsteps);
  CProxy_Node(CpvAccess(BOCclass_group).node).scriptParam(msg);
  Node::Object()->enableScriptBarrier();
  script->suspend();

  script->scriptBarrier.publish(script->barrierStep++,1);
  script->suspend();

  msg = new ScriptParamMsg;
  sprintf(msg->param,"firsttimestep");
  sprintf(msg->value,"%d",simParams->N);
  CProxy_Node(CpvAccess(BOCclass_group).node).scriptParam(msg);
  Node::Object()->enableScriptBarrier();
  script->suspend();

  return TCL_OK;
}

int ScriptTcl::Tcl_move(ClientData clientData,
	Tcl_Interp *interp, int argc, char *argv[]) {
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
  ScriptTcl *script = (ScriptTcl *)clientData;

  iout << "TCL: Moving atom " << atomid << " ";
  if ( moveto ) iout << "to"; else iout << "by";
  iout << " " << Vector(x,y,z) << ".\n" << endi;

  MoveAtomMsg *msg = new MoveAtomMsg;
  msg->atomid = atomid - 1;
  msg->moveto = moveto;
  msg->coord = Vector(x,y,z);
  CProxy_PatchMgr(CpvAccess(BOCclass_group).patchMgr).moveAtom(msg);
  Node::Object()->enableScriptBarrier();
  script->suspend();

  return TCL_OK;
}

int ScriptTcl::Tcl_output(ClientData clientData,
        Tcl_Interp *interp, int argc, char *argv[]) {

  iout << "TCL: Triggering file output.\n" << endi;

  ScriptTcl *script = (ScriptTcl *)clientData;
  script->scriptBarrier.publish(script->barrierStep++,2);
  script->suspend();

  return TCL_OK;
}

#endif  // NAMD_TCL


ScriptTcl::ScriptTcl() : scriptBarrier(scriptBarrierTag) {
  DebugM(3,"Constructing ScriptTcl\n");
#ifdef NAMD_TCL
  interp = 0;
#endif
}

ScriptTcl::~ScriptTcl() {
  DebugM(3,"Destructing ScriptTcl\n");
#ifdef NAMD_TCL
  if ( interp ) Tcl_DeleteInterp(interp);
#endif
}

// Invoked by thread
void ScriptTcl::threadRun(ScriptTcl* arg)
{
    arg->algorithm();
}

// Invoked by Node::run()
void ScriptTcl::run()
{
    // create a Thread and invoke it
    DebugM(4, "::run() - this = " << this << "\n" );
    thread = CthCreate((CthVoidFn)&(threadRun),(void*)(this),TCL_STK_SZ);
    CthSetStrategyDefault(thread);
    // awaken(); // triggered by sequencers
}

void ScriptTcl::algorithm() {
  DebugM(4,"Running ScriptTcl\n");

  barrierStep = 0;

#ifdef NAMD_TCL
  // Create interpreter
  interp = Tcl_CreateInterp();
  if (Tcl_Init(interp) == TCL_ERROR) {
    CkPrintf("Tcl startup error: %\n", interp->result);
  }
  Tcl_CreateCommand(interp, "print", Tcl_print,
    (ClientData) NULL, (Tcl_CmdDeleteProc *) NULL);
  Tcl_CreateCommand(interp, "param", Tcl_param,
    (ClientData) this, (Tcl_CmdDeleteProc *) NULL);
  Tcl_CreateCommand(interp, "run", Tcl_run,
    (ClientData) this, (Tcl_CmdDeleteProc *) NULL);
  Tcl_CreateCommand(interp, "move", Tcl_move,
    (ClientData) this, (Tcl_CmdDeleteProc *) NULL);
  Tcl_CreateCommand(interp, "output", Tcl_output,
    (ClientData) this, (Tcl_CmdDeleteProc *) NULL);

  // Get the script
  StringList *script = Node::Object()->configList->find("tclScript");

  for ( ; script; script = script->next ) {
    int code;
    if ( script->data[0] == '{' ) code = Tcl_Eval(interp,script->data+1);
    else code = Tcl_EvalFile(interp,script->data);
    if (*interp->result != 0) CkPrintf("TCL: %s\n",interp->result);
    if (code != TCL_OK) NAMD_die("TCL error in script!");
  }

#else

  NAMD_die("Sorry, TCL scripting not available; built without TCL.");

#endif

  scriptBarrier.publish(barrierStep++,0);  // terminate sequencers
  suspend();

}


/***************************************************************************
 * RCS INFORMATION:
 *
 *	$RCSfile $
 *	$Author $	$Locker:  $		$State: Exp $
 *	$Revision: 1.4 $	$Date: 1999/08/30 14:53:32 $
 *
 ***************************************************************************
 * REVISION HISTORY:
 *
 * $Log: ScriptTcl.C,v $
 * Revision 1.4  1999/08/30 14:53:32  jim
 * Added error messages for TCL features if TCL unavailable.
 *
 * Revision 1.3  1999/08/11 16:53:10  jim
 * Added move command to TCL scripting.
 *
 * Revision 1.2  1999/06/21 16:15:36  jim
 * Improved scripting, run now ends and generates output.
 *
 * Revision 1.1  1999/05/27 18:38:58  jim
 * Files to implement general Tcl scripting.
 *
 *
 *
 ***************************************************************************/
