/**
***  Copyright (c) 1995, 1996, 1997, 1998, 1999, 2000 by
***  The Board of Trustees of the University of Illinois.
***  All rights reserved.
**/

#include "Measure.h"
#include "Node.h"
#include "Parameters.h"
#include "Molecule.h"

#ifdef NAMD_TCL

int Measure::wrapCommand(ClientData clientData,
        Tcl_Interp *interp, int argc, char *argv[]) {

/*
  double (wraped *)(Vector*, Molecule*, Parameters*);
  wraped = 
  Node *node = Node::Object();
  Molecule *molecule = node->molecule;
  Parameters *parameters = node->parameters;
  Vector *coordinates = node->coords;

  double result = *wraped(coordinates,molecule,parameters);

*/
  return TCL_OK;
}

int Tcl_centerOfNumber(ClientData, Tcl_Interp *interp, int argc, char *argv[]) {

  Node *node = Node::Object();
  Molecule *molecule = node->molecule;
  Parameters *parameters = node->parameters;
  Vector *coordinates = node->coords;
  int numAtoms = molecule->numAtoms;

  int number = 0;
  Vector center;
  for( int i = 0; i < numAtoms; ++i ) {
    number += 1;
    center += coordinates[i];
  }
  center /= number;

  char s[1024];
  sprintf(s,"{ %g %g %g }", center.x, center.y, center.z);
  Tcl_SetResult(interp,s,TCL_VOLATILE);

  return TCL_OK;
}

int Tcl_centerOfMass(ClientData, Tcl_Interp *interp, int argc, char *argv[]) {

  Node *node = Node::Object();
  Molecule *molecule = node->molecule;
  Parameters *parameters = node->parameters;
  Vector *coordinates = node->coords;
  int numAtoms = molecule->numAtoms;

  Vector center;
  BigReal totalMass = 0;
  for( int i = 0; i < numAtoms; ++i ) {
    BigReal mass = molecule->atommass(i);
    totalMass += mass;
    center += mass * coordinates[i];
  }
  center /= totalMass;

  char s[1024];
  sprintf(s,"{ %g %g %g }", center.x, center.y, center.z);
  Tcl_SetResult(interp,s,TCL_VOLATILE);

  return TCL_OK;
}

void Measure::createCommands(Tcl_Interp *interp) {
  Tcl_CreateCommand(interp, "centerOfNumber", Tcl_centerOfNumber,
    (ClientData) NULL, (Tcl_CmdDeleteProc *) NULL);
  Tcl_CreateCommand(interp, "centerOfMass", Tcl_centerOfMass,
    (ClientData) NULL, (Tcl_CmdDeleteProc *) NULL);
}

void Measure::deleteCommands(Tcl_Interp *interp) {
  Tcl_DeleteCommand(interp, "centerOfNumber");
  Tcl_DeleteCommand(interp, "centerOfMass");
}

#endif

