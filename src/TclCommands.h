/***************************************************************************
 *cr                                                                       
 *cr            (C) Copyright 1995 The Board of Trustees of the           
 *cr                        University of Illinois                       
 *cr                         All Rights Reserved                        
 *cr                                                                   
 ***************************************************************************/

/***************************************************************************
 * DESCRIPTION:
 *  prototypes for VMD<->Tcl functions
 *
 ***************************************************************************/
#ifndef TCLCOMMANDS_H
#define TCLCOMMANDS_H

#ifdef NAMD_TCL

#include <tcl.h>

//forward definition
class Matrix4;
class AtomSel;

// override some vector functions for speed
int proc_vecadd(ClientData, Tcl_Interp *interp, int argc, char *argv[]);
int proc_vecsub(ClientData, Tcl_Interp *interp, int argc, char *argv[]);
int proc_vecscale(ClientData, Tcl_Interp *interp, int argc, char *argv[]);
int proc_transoffset(ClientData, Tcl_Interp *interp, int argc, char *argv[]);
int proc_transmult(ClientData, Tcl_Interp *interp, int argc, char *argv[]);
int proc_vectrans(ClientData, Tcl_Interp *interp, int argc, char *argv[]);

// get a matrix from a string; 
// returns TCL_OK if good
// If bad, returns TCL_ERROR and sets the interp->result to the error message
// The name of the function should be passed in 'fctn' so the error message
// can be constructed correctly
int tcl_get_matrix(char *fctn, Tcl_Interp *interp,
		   char *s, Matrix4 *mat);

// append the matrix information to the interp->result field
void tcl_append_matrix(Tcl_Interp *interp, const Matrix4 &mat);

// get a vector -- YOU must delete the vector, if successful
// returns TCL_OK if good
// If bad, returns TCL_ERROR and sets the interp->result to the error message
// The name of the function should be passed in 'fctn' so the error message
// can be constructed correctly
int tcl_get_vector(char *fctn, Tcl_Interp *interp, 
			  char *s, int *num, float **result);

// append the double to the end of the string s
// return the end of the string
// optionally add a space at the beginning
char *tcl_append_double(char *s, double f, int add_space = 0);

#endif
#endif
