/**
***  Copyright (c) 1995, 1996, 1997, 1998, 1999, 2000 by
***  The Board of Trustees of the University of Illinois.
***  All rights reserved.
**/

/*
   Forwards atoms to master node for force evaluation.
*/

#include "Node.h"
#include "ComputeTcl.h"
#include "ComputeGlobal.h"
#include "ComputeGlobalMsgs.h"
#include "Molecule.h"
#include "ReductionMgr.h"
#include "ComputeMgr.h"
#include "ComputeMgr.decl.h"
#include "ScriptTcl.h"
#include <stdio.h>

#ifdef NAMD_TCL
#include <tcl.h>
#endif
#include "TclCommands.h"

// #define DEBUGM
#define MIN_DEBUG_LEVEL 4
#include "Debug.h"


#ifdef NAMD_TCL
int ComputeTcl::Tcl_print(ClientData,
	Tcl_Interp *, int argc, char *argv[]) {
  int arglen = 1;  int ai;
  for (ai=1; ai<argc; ++ai) { arglen += strlen(argv[ai]) + 1; }
  char *buf = new char[arglen];  *buf = 0;
  for (ai=1; ai<argc; ++ai) { strcat(buf,argv[ai]); strcat(buf," "); }
  ai = strlen(buf);  if ( ai ) buf[ai-1] = 0;
  CkPrintf("TCL: %s\n",buf);
  delete [] buf;
  return TCL_OK;
}


int ComputeTcl::Tcl_atomid(ClientData clientData,
	Tcl_Interp *interp, int argc, char *argv[]) {
  if (argc != 4) {
    Tcl_SetResult(interp,"wrong # args",TCL_VOLATILE);
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
    Tcl_SetResult(interp,"atom not found",TCL_VOLATILE);
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
    Tcl_SetResult(interp,"wrong # args",TCL_VOLATILE);
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


int ComputeTcl::Tcl_addgroup(ClientData clientData,
	Tcl_Interp *interp, int argc, char *argv[]) {
  if (argc != 2) {
    Tcl_SetResult(interp,"wrong # args",TCL_VOLATILE);
    return TCL_ERROR;
  }

  AtomIDList *gdef = (AtomIDList *)clientData;
  AtomIDList::iterator g_i, g_e;
  g_i = gdef->begin();  g_e = gdef->end();
  int gcount = 1;
  for ( ; g_i != g_e; ++g_i ) if ( *g_i == -1 ) ++gcount;

  int listc, i;  char **listv;

  if (Tcl_SplitList(interp,argv[1],&listc,&listv) != TCL_OK) {
    return TCL_ERROR;
  }
  for ( i = 0; i < listc; ++i ) {
    int atomid;
    if (Tcl_GetInt(interp,listv[i],&atomid) != TCL_OK) {
      gdef->add(-1);  // sentinel - this should use a delete instead
      Tcl_Free((char*) listv);
      return TCL_ERROR;
    }
    gdef->add(atomid-1);
  }
  gdef->add(-1);  // sentinel
  Tcl_Free((char*) listv);

  char s[10];  sprintf(s,"g%d",gcount);
  Tcl_SetResult(interp,s,TCL_VOLATILE);

  DebugM(4,"Group " << s << " added to config list\n");
  return TCL_OK;
}


int ComputeTcl::Tcl_reconfig(ClientData clientData,
	Tcl_Interp *interp, int argc, char **) {
  if (argc != 1) {
    Tcl_SetResult(interp,"wrong # args",TCL_VOLATILE);
    return TCL_ERROR;
  }
  int *reconfig = (int *)clientData;
  *reconfig = 1;
  DebugM(4,"Reconfiguration turned on\n");
  return TCL_OK;
}


int ComputeTcl::Tcl_loadcoords(ClientData clientData,
	Tcl_Interp *interp, int objc, Tcl_Obj * const objv[]) {
  if (objc != 2) {
    Tcl_SetResult(interp,"wrong # args",TCL_VOLATILE);
    return TCL_ERROR;
  }
  Tcl_Obj * const vname = objv[1];
  ComputeTcl *self = (ComputeTcl *)clientData;
  AtomIDList::iterator a_i = self->aid.begin();
  AtomIDList::iterator a_e = self->aid.end();
  PositionList::iterator p_i = self->p.begin();
  for ( ; a_i != a_e; ++a_i, ++p_i ) {
    Tcl_Obj *newlist = Tcl_NewListObj(0, NULL);
    Tcl_Obj *arrkey = Tcl_NewIntObj((int)((*a_i)+1));
    
    Tcl_ListObjAppendElement(interp, newlist, 
      Tcl_NewDoubleObj((double)((*p_i).x)));
    Tcl_ListObjAppendElement(interp, newlist, 
      Tcl_NewDoubleObj((double)((*p_i).y)));
    Tcl_ListObjAppendElement(interp, newlist, 
      Tcl_NewDoubleObj((double)((*p_i).z)));
   
    if (!Tcl_ObjSetVar2(interp, vname, arrkey, newlist, 0)) {
      NAMD_die("TCL error in global force calculation!");
      return TCL_ERROR;
    }
  }
  PositionList::iterator c_i = self->gcom.begin();
  PositionList::iterator c_e = self->gcom.end();
  int gcount = 1;
  for ( ; c_i != c_e; ++c_i, ++gcount ) {
    Tcl_Obj *newlist = Tcl_NewListObj(0, NULL);
    char buf[10];
    sprintf(buf, "g%d", gcount);
    Tcl_Obj *arrkey = Tcl_NewStringObj(buf, -1);
 
    Tcl_ListObjAppendElement(interp, newlist,
      Tcl_NewDoubleObj((double)((*c_i).x)));
    Tcl_ListObjAppendElement(interp, newlist,
      Tcl_NewDoubleObj((double)((*c_i).y)));
    Tcl_ListObjAppendElement(interp, newlist,
      Tcl_NewDoubleObj((double)((*c_i).z)));
   
    if (!Tcl_ObjSetVar2(interp, vname, arrkey, newlist, 0)) {
      NAMD_die("TCL error in global force calculation!");
      return TCL_ERROR;
    }
  }
  return TCL_OK;
}


int ComputeTcl::Tcl_loadmasses(ClientData clientData,
	Tcl_Interp *interp, int objc, Tcl_Obj * const objv[]) {
  if (objc != 2) {
    Tcl_SetResult(interp,"wrong # args",TCL_VOLATILE);
    return TCL_ERROR;
  }
  Tcl_Obj * const vname = objv[1];
  ComputeTcl *self = (ComputeTcl *)clientData;
  Molecule *mol = Node::Object()->molecule;
  AtomIDList::iterator a_i = self->aid.begin();
  AtomIDList::iterator a_e = self->aid.end();
  for ( ; a_i != a_e; ++a_i) {
    if (!Tcl_ObjSetVar2(interp, vname,
                        Tcl_NewIntObj((int)((*a_i)+1)),
                        Tcl_NewDoubleObj((double)(mol->atommass(*a_i))),
                        0)) {
      NAMD_die("TCL error in global force calculation!");
      return TCL_ERROR;
    }
  }
  ResizeArray<BigReal>::iterator g_i, g_e;
  g_i = self->gmass.begin();  g_e = self->gmass.end();
  int gcount = 1;
  for ( ; g_i != g_e; ++g_i, ++gcount) {
    char buf[10];
    sprintf(buf, "g%d", gcount);
    if (!Tcl_ObjSetVar2(interp, vname,
                        Tcl_NewStringObj(buf, -1),
                        Tcl_NewDoubleObj((double)(*g_i)),
                        0)) {
      NAMD_die("TCL error in global force calculation!");
      return TCL_ERROR;
    }
  }
  return TCL_OK;
}


int ComputeTcl::Tcl_addforce(ClientData clientData,
	Tcl_Interp *interp, int objc, Tcl_Obj * const objv[]) {
  if (objc != 3) {
    Tcl_SetResult(interp,"wrong # args",TCL_VOLATILE);
    return TCL_ERROR;
  }
  Tcl_Obj **force;  int fnum;  int atomid;  double x, y, z;
  int isgroup = 0;
  char *id = Tcl_GetStringFromObj(objv[1], NULL); 
  if ( id[0] == 'g' ) {
    isgroup = 1;
    if ( Tcl_GetInt(interp,id+1,&atomid) != TCL_OK ) return TCL_ERROR;
  } else {
    if ( Tcl_GetInt(interp,id,&atomid) != TCL_OK ) return TCL_ERROR;
  }
  if (Tcl_ListObjGetElements(interp, objv[2], &fnum, &force) != TCL_OK) {
    return TCL_ERROR;
  }
  if ( (fnum != 3) ||
       (Tcl_GetDoubleFromObj(interp, force[0],&x) != TCL_OK) ||
       (Tcl_GetDoubleFromObj(interp, force[1],&y) != TCL_OK) ||
       (Tcl_GetDoubleFromObj(interp, force[2],&z) != TCL_OK) ) {
    Tcl_SetResult(interp,"force not a vector",TCL_VOLATILE);
    return TCL_ERROR;
  }
  ComputeGlobalResultsMsg *msg = (ComputeGlobalResultsMsg *)clientData;
  if ( isgroup ) {
    msg->gforce.item(atomid-1) += Vector(x,y,z);
  } else {
    msg->aid.add(atomid-1);
    msg->f.add(Vector(x,y,z));
  }
  DebugM(4,"Atom ID " << atomid << " added to force list\n");
  return TCL_OK;
}
#endif


ComputeTcl::ComputeTcl(ComputeMgr *c) : ComputeGlobalMaster(c) {
  DebugM(3,"Constructing ComputeTcl\n");
#ifdef NAMD_TCL
  interp = 0;
#endif
}

ComputeTcl::~ComputeTcl() {
  DebugM(3,"Destructing ComputeTcl\n");
#ifdef NAMD_TCL
/*
  if ( interp ) Tcl_DeleteInterp(interp);
*/
#endif
}


void ComputeTcl::initialize() {
  DebugM(4,"Initializing master\n");

  ComputeGlobalConfigMsg *msg = new ComputeGlobalConfigMsg;
  msg->tag = tag;
   

#ifdef NAMD_TCL
  interp = Node::Object()->getScript()->interp;
  Tcl_CreateCommand(interp, "atomid", Tcl_atomid,
    (ClientData) (Node::Object()->molecule), (Tcl_CmdDeleteProc *) NULL);
  Tcl_CreateCommand(interp, "vecadd", proc_vecadd,
    (ClientData) NULL, (Tcl_CmdDeleteProc *) NULL);
  Tcl_CreateCommand(interp, "vecsub", proc_vecsub,
    (ClientData) NULL, (Tcl_CmdDeleteProc *) NULL);
  Tcl_CreateCommand(interp, "vecscale", proc_vecscale,
    (ClientData) NULL, (Tcl_CmdDeleteProc *) NULL);

  // Call interpreter to determine requested atoms
  Tcl_CreateCommand(interp, "addatom", Tcl_addatom,
    (ClientData) &(msg->aid), (Tcl_CmdDeleteProc *) NULL);
  Tcl_CreateCommand(interp, "addgroup", Tcl_addgroup,
    (ClientData) &(msg->gdef), (Tcl_CmdDeleteProc *) NULL);

  // Get the script
  StringList *script = Node::Object()->configList->find("tclForcesScript");

  for ( ; script; script = script->next ) {
    int code;
    if ( strstr(script->data,"\n") ) {
       code = Tcl_Eval(interp,script->data);
    }
    else code = Tcl_EvalFile(interp,script->data);
    char *result = Tcl_GetStringResult(interp);
    if (*result != 0) CkPrintf("TCL: %s\n",result);
    if (code != TCL_OK) {
      char *errorInfo = Tcl_GetVar(interp,"errorInfo",0);
      NAMD_die(errorInfo);
    }
  }

  Tcl_DeleteCommand(interp, "addatom");
  Tcl_DeleteCommand(interp, "addgroup");
#else

  NAMD_die("Sorry, tclForces is not available; built without TCL.");

#endif

  storedefs(msg->gdef);

  // Send config to clients
  comm->sendComputeGlobalConfig(msg);
}


void ComputeTcl::calculate() {
  DebugM(4,"Calculating forces on master\n");

  ComputeGlobalResultsMsg *msg = new ComputeGlobalResultsMsg;
  msg->tag = tag;
  msg->gforce.resize(gmass.size());
  msg->gforce.setall(Vector(0,0,0));

#ifdef NAMD_TCL
  // Call interpreter to calculate forces
  Tcl_CreateObjCommand(interp, (char *)"loadcoords", Tcl_loadcoords,
    (ClientData) this, (Tcl_CmdDeleteProc *) NULL);
  Tcl_CreateObjCommand(interp, (char *)"loadmasses", Tcl_loadmasses,
    (ClientData) this, (Tcl_CmdDeleteProc *) NULL);
  Tcl_CreateObjCommand(interp, (char *)"addforce", Tcl_addforce,
    (ClientData) msg, (Tcl_CmdDeleteProc *) NULL);
  Tcl_CreateCommand(interp, (char *)"reconfig", Tcl_reconfig,
    (ClientData) &(msg->reconfig), (Tcl_CmdDeleteProc *) NULL);
  Tcl_CreateCommand(interp, (char *)"addatom", Tcl_addatom,
    (ClientData) &(msg->newaid), (Tcl_CmdDeleteProc *) NULL);
  Tcl_CreateCommand(interp, (char *)"addgroup", Tcl_addgroup,
    (ClientData) &(msg->newgdef), (Tcl_CmdDeleteProc *) NULL);

  char cmd[129];  int code;
  strcpy(cmd,"calcforces");  code = Tcl_Eval(interp,cmd);
  char *result = Tcl_GetStringResult(interp);
  if (*result != 0) CkPrintf("TCL: %s\n",result);
  if (code != TCL_OK) {
    char *errorInfo = Tcl_GetVar(interp,"errorInfo",0);
    NAMD_die(errorInfo);
  }

  Tcl_DeleteCommand(interp, "loadcoords");
  Tcl_DeleteCommand(interp, "loadmasses");
  Tcl_DeleteCommand(interp, "addforce");
  Tcl_DeleteCommand(interp, "reconfig");
  Tcl_DeleteCommand(interp, "addatom");
  Tcl_DeleteCommand(interp, "addgroup");
#endif

  // Send results to clients
  DebugM(3,"Sending results (" << msg->aid.size() << " forces) on master\n");
  if ( msg->reconfig ) {
    storedefs(msg->newgdef);
    DebugM(4,"Sending new configuration (" <<
			msg->newaid.size() << " atoms) on master\n");
  }
  comm->sendComputeGlobalResults(msg);
}

