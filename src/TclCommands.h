/**
***  Copyright (c) 1995, 1996, 1997, 1998, 1999, 2000 by
***  The Board of Trustees of the University of Illinois.
***  All rights reserved.
**/

#ifndef TCLCOMMANDS_H
#define TCLCOMMANDS_H

#ifdef NAMD_TCL

#include <tcl.h>

int proc_vecadd(ClientData, Tcl_Interp *interp, int argc, char *argv[]);
int proc_vecsub(ClientData, Tcl_Interp *interp, int argc, char *argv[]);
int proc_vecscale(ClientData, Tcl_Interp *interp, int argc, char *argv[]);

int tcl_get_vector(char *fctn, Tcl_Interp *interp, 
			  char *s, int *num, float **result);

#endif
#endif
