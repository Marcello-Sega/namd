/**
***  Copyright (c) 1995, 1996, 1997, 1998, 1999, 2000 by
***  The Board of Trustees of the University of Illinois.
***  All rights reserved.
**/

#include <stdlib.h> 
#include <malloc.h>
#include <errno.h>
#include "TclCommands.h"

#ifdef NAMD_TCL

#include <tcl.h>

#define SIMPLE_TCL_OPT(string,result)       \
if (!strcmp(argv[1], string)) {             \
  Tcl_AppendResult(interp, result, NULL);   \
  return TCL_OK;                            \
}


/***************** override some of the vector routines for speed ******/
/* These should be the exact C equivalent to the corresponding Tcl    */
/* vector commands */

// Function:  vecadd v1 v2 {v3 ...}
//  Returns: the sum of vectors; must all be the same length
//  The increase in speed from Tcl to C++ is 4561 / 255 == 18 fold
int proc_vecadd(ClientData, Tcl_Interp *interp, int argc, 
		       char *argv[])
{
  if (argc == 1) {
    Tcl_SetResult(interp,"no value given for parameter \"x\" to \"vecadd\"",TCL_VOLATILE);
    return TCL_ERROR;
  }
  if (argc == 2) {
    Tcl_SetResult(interp,"no value given for parameter \"y\" to \"vecadd\"",TCL_VOLATILE);
    return TCL_ERROR;
  }
  int num;
  char **data;
  if (Tcl_SplitList(interp, argv[1], &num, &data) != TCL_OK) {
    return TCL_ERROR;
  }
  double *sum = new double[num];
  int i;
  for (i=0; i<num; i++) {
    if (Tcl_GetDouble(interp, data[i], sum+i) != TCL_OK) {
      delete [] sum;
      Tcl_Free((char*) data);
      return TCL_ERROR;
    }
  }
  Tcl_Free((char*) data);
  // do the sums on the rest
  int num2;
  for (int term=2; term < argc; term++) {
    if (Tcl_SplitList(interp, argv[term], &num2, &data) != TCL_OK) {
      delete [] sum;
      return TCL_ERROR;
    }
    if (num != num2) {
      Tcl_SetResult(interp,"vecadd: two vectors don't have the same size",TCL_VOLATILE);
      delete [] sum;
      Tcl_Free((char*) data);
      return TCL_ERROR;
    }
    for (i=0; i<num; i++) {
      double df;
      if (Tcl_GetDouble(interp, data[i], &df) != TCL_OK) {
	delete [] sum;
	Tcl_Free((char*) data);
	return TCL_ERROR;
      }
      sum[i] += df;
    }
  }

  // and return the result
  char s[TCL_DOUBLE_SPACE];
  for (i=0; i<num; i++) {
    Tcl_PrintDouble(interp, sum[i], s);
    Tcl_AppendElement(interp, s);
  }
  Tcl_Free((char*) data);
  delete [] sum;
  return TCL_OK;
}

// Function:  vecsub  v1 v2
//  Returns:   v1 - v2
int proc_vecsub(ClientData, Tcl_Interp *interp, int argc, char *argv[])
{
  if (argc == 1) {
    Tcl_SetResult(interp,"no value given for parameter \"x\" to \"vecsub\"",TCL_VOLATILE);
    return TCL_ERROR;
  }
  if (argc == 2) {
    Tcl_SetResult(interp,"no value given for parameter \"y\" to \"vecsub\"",TCL_VOLATILE);
    return TCL_ERROR;
  }
  int num1, num2;
  float *data1, *data2;
  if (tcl_get_vector("vecsub: ", interp, argv[1], &num1, &data1) != TCL_OK) {
    return TCL_ERROR;
  }
  if (tcl_get_vector("vecsub: ", interp, argv[2], &num2, &data2) != TCL_OK) {
    delete [] data1;
    return TCL_ERROR;
  }
  if (num1 != num2) {
    Tcl_SetResult(interp,"vecadd: two vectors don't have the same size",TCL_VOLATILE);
    delete [] data1;
    delete [] data2;
    return TCL_ERROR;
  }
  // do the subtraction and return the result
  char s[TCL_DOUBLE_SPACE];
  for (int i=0; i<num1; i++) {
    Tcl_PrintDouble(interp, data1[i] - data2[i], s);
    Tcl_AppendElement(interp, s);
  }
  delete [] data1;
  delete [] data2;
  return TCL_OK;
}


// Function: vecscale
//  Returns: scalar * vector or vector * scalar
// speedup is 1228/225 = 5.5 fold
int proc_vecscale(ClientData, Tcl_Interp *interp, int argc, 
		       char *argv[])
{
  if (argc == 1) {
    Tcl_SetResult(interp,"no value given for parameter \"c\" to \"vecscale\"",TCL_VOLATILE);
    return TCL_ERROR;
  }
    
  if (argc == 2) {
    Tcl_SetResult(interp,"no value given for parameter \"v\" to \"vecscale\"",TCL_VOLATILE);
    return TCL_ERROR;
  }
  if (argc != 3) {
    Tcl_SetResult(interp,"called \"vecscale\" with too many arguments",TCL_VOLATILE);
    return TCL_ERROR;
  }
    
  int num1, num2;
  char **data1, **data2;
  if (Tcl_SplitList(interp, argv[1], &num1, &data1) != TCL_OK) {
    return TCL_ERROR;
  }
  if (Tcl_SplitList(interp, argv[2], &num2, &data2) != TCL_OK) {
    Tcl_Free((char*) data1);
    return TCL_ERROR;
  }
  int result = TCL_OK;
  if (num1 == 0 || num2 == 0) {
    result = TCL_ERROR;
    Tcl_SetResult(interp,"vecscale: parameters must have data",TCL_VOLATILE);
  } else if (num1 != 1 && num2 != 1) {
    result = TCL_ERROR;
    Tcl_SetResult(interp,"vecscale: one parameter must be a scalar value",TCL_VOLATILE);
  } else {
    char *scalar, **vector;
    int num;
    if (num1 == 1) {
      scalar = data1[0];
      vector = data2;
      num = num2;
    } else {
      scalar = data2[0];
      vector = data1;
      num = num1;
    }
    char s[TCL_DOUBLE_SPACE];
    double val1, val2;
    if (Tcl_GetDouble(interp, scalar, &val1) != TCL_OK) {
      result = TCL_ERROR;
    } else {
      for (int i=0; i<num; i++) {
	if (Tcl_GetDouble(interp, vector[i], &val2) != TCL_OK) {
	  Tcl_SetResult(interp,"vecscale: vector contains a non-number",TCL_VOLATILE);
	  result = TCL_ERROR;
	  break;
	}
	Tcl_PrintDouble(interp, val1 * val2, s);
	Tcl_AppendElement(interp, s);
      }
    }
  }
  Tcl_Free((char*) data1);
  Tcl_Free((char*) data2);
  return result;
}


// Given a string with a vector in it, get the vector
// YOU must delete [] the vector (in "result") when finished
// returns TCL_OK if good
// If bad, returns TCL_ERROR and sets the result to the error message
// The name of the function should be passed in 'fctn' so the error message
// can be constructed correctly
int tcl_get_vector(char *fctn, Tcl_Interp *interp, 
			  char *s, int *num, float **result)
{
  *result = NULL;
  *num = 0;
  char **data;
  if (Tcl_SplitList(interp, s, num, &data) != TCL_OK) { // is a list
    Tcl_AppendResult(interp, fctn, ": badly formed vector", NULL);
    return TCL_ERROR;
  }
  *result = new float[*num];
  int ret_val = TCL_OK;
  double tmp;
  for (int i=0; i<*num; i++) {
    if (Tcl_GetDouble(interp, data[i], &tmp) != TCL_OK) {  // of numbers
      Tcl_SetResult(interp,"non-numeric in vector",TCL_VOLATILE);
      ret_val = TCL_ERROR;
    } else {
      (*result)[i] = tmp;
    }
  }
  Tcl_Free((char*) data);
  if (ret_val == TCL_ERROR) {
    delete [] (*result);
    *result = NULL;
  }
  return ret_val;
}

#endif
