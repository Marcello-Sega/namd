/***************************************************************************
 *cr                                                                       
 *cr            (C) Copyright 1995 The Board of Trustees of the           
 *cr                        University of Illinois                       
 *cr                         All Rights Reserved                        
 *cr                                                                   
 ***************************************************************************/

/***************************************************************************
 * RCS INFORMATION:
 *
 *	$RCSfile: TclCommands.C,v $
 *	$Author: jim $	$Locker:  $		$State: Exp $
 *	$Revision: 1.5 $	$Date: 1999/06/28 19:49:57 $
 *
 ***************************************************************************
 * DESCRIPTION:
 *   Tcl <--> VMD interface commands used for the analysis and 
 * manipulation of structures
 *
 ***************************************************************************/

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
    interp -> result = "no value given for parameter \"x\" to \"vecadd\"";
    return TCL_ERROR;
  }
  if (argc == 2) {
    interp -> result = "no value given for parameter \"y\" to \"vecadd\"";
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
      free(data);
      return TCL_ERROR;
    }
  }
  free(data);
  // do the sums on the rest
  int num2;
  for (int term=2; term < argc; term++) {
    if (Tcl_SplitList(interp, argv[term], &num2, &data) != TCL_OK) {
      delete [] sum;
      return TCL_ERROR;
    }
    if (num != num2) {
      interp -> result = "vecadd: two vectors don't have the same size";
      delete [] sum;
      free(data);
      return TCL_ERROR;
    }
    for (i=0; i<num; i++) {
      double df;
      if (Tcl_GetDouble(interp, data[i], &df) != TCL_OK) {
	delete [] sum;
	free(data);
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
  free(data);
  delete [] sum;
  return TCL_OK;
}

// Function:  vecsub  v1 v2
//  Returns:   v1 - v2
int proc_vecsub(ClientData, Tcl_Interp *interp, int argc, char *argv[])
{
  if (argc == 1) {
    interp -> result = "no value given for parameter \"x\" to \"vecsub\"";
    return TCL_ERROR;
  }
  if (argc == 2) {
    interp -> result = "no value given for parameter \"y\" to \"vecsub\"";
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
    interp -> result = "vecadd: two vectors don't have the same size";
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
  return TCL_OK;
}


// Function: vecscale
//  Returns: scalar * vector or vector * scalar
// speedup is 1228/225 = 5.5 fold
int proc_vecscale(ClientData, Tcl_Interp *interp, int argc, 
		       char *argv[])
{
  if (argc == 1) {
    interp -> result = "no value given for parameter \"c\" to \"vecscale\"";
    return TCL_ERROR;
  }
    
  if (argc == 2) {
    interp -> result = "no value given for parameter \"v\" to \"vecscale\"";
    return TCL_ERROR;
  }
  if (argc != 3) {
    interp -> result = "called \"vecscale\" with too many arguments";
    return TCL_ERROR;
  }
    
  int num1, num2;
  char **data1, **data2;
  if (Tcl_SplitList(interp, argv[1], &num1, &data1) != TCL_OK) {
    return TCL_ERROR;
  }
  if (Tcl_SplitList(interp, argv[2], &num2, &data2) != TCL_OK) {
    free(data1);
    return TCL_ERROR;
  }
  int result = TCL_OK;
  if (num1 == 0 || num2 == 0) {
    result = TCL_ERROR;
    interp -> result = "vecscale: parameters must have data";
  } else if (num1 != 1 && num2 != 1) {
    result = TCL_ERROR;
    interp -> result = "vecscale: one parameter must be a scalar value";
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
	  interp -> result = "vecscale: vector contains a non-number";
	  result = TCL_ERROR;
	  break;
	}
	Tcl_PrintDouble(interp, val1 * val2, s);
	Tcl_AppendElement(interp, s);
      }
    }
  }
  free(data1);
  free(data2);
  return result;
}

/*

//  Function: transoffset
//   Returns: the transformation correspoding to a vector offset
int proc_transoffset(ClientData, Tcl_Interp *interp, int argc, 
		     char *argv[])
{
  if (argc != 2) {
    interp -> result = "transoffset: takes one parameter, an offset vector";
    return TCL_ERROR;
  }
  // get the vector
  int num;
  float *data;
  if (tcl_get_vector("transoffset: ", interp, argv[1], &num, &data) != 
      TCL_OK) {
    return TCL_ERROR;
  }
  // don't check the size (the script version doesn't)
  Matrix4 t;
  switch (num) {
  case 3: t.mat[3][2] = data[2];  // YES, these fall through!
  case 2: t.mat[3][1] = data[1];
  case 1: t.mat[3][0] = data[0];
  case 0: break;
  }
  delete [] data;
  tcl_append_matrix(interp, t);
  return TCL_OK;
  
}

/// Given a string with a matrix in it, return the matrix
// returns TCL_OK if good
// If bad, returns TCL_ERROR and sets the interp->result to the error message
// The name of the function should be passed in 'fctn' so the error message
// can be constructed correctly
int tcl_get_matrix(char *fctn, Tcl_Interp *interp, 
			  char *s, Matrix4 *mat)
{ 
  int num_rows;
  char **data_rows;
  if (Tcl_SplitList(interp, s, &num_rows, &data_rows) != TCL_OK) {
    sprintf(interp -> result, "%s: badly formed matrix", fctn);
    return TCL_ERROR;
  }
  if (num_rows != 4) {
    free(data_rows);
    sprintf(interp -> result, "%s: need a 4x4 matrix", fctn);
    return TCL_ERROR;
  }
  int num_row[4];
  char **data_row[4];
  data_row[0] = data_row[1] = data_row[2] = data_row[3] = NULL;
  if (Tcl_SplitList(interp, data_rows[0], num_row+0, data_row+0) != TCL_OK ||
      num_row[0] != 4 ||
      Tcl_SplitList(interp, data_rows[1], num_row+1, data_row+1) != TCL_OK ||
      num_row[1] != 4 ||
      Tcl_SplitList(interp, data_rows[2], num_row+2, data_row+2) != TCL_OK ||
      num_row[2] != 4 ||
      Tcl_SplitList(interp, data_rows[3], num_row+3, data_row+3) != TCL_OK ||
      num_row[3] != 4) {
    free(data_rows);
    if (data_row[0]) free(data_row[0]);
    if (data_row[1]) free(data_row[1]);
    if (data_row[2]) free(data_row[2]);
    if (data_row[3]) free(data_row[3]);
    Tcl_AppendResult(interp, fctn, ": poorly formed matrix", NULL);
    return TCL_ERROR;
  }
  free(data_rows);
  // now get the numbers
  double tmp;
  int ret_val = TCL_OK;
  for (int i=0; i<4; i++) {
    for (int j=0; j<4; j++) {
      if (Tcl_GetDouble(interp, data_row[i][j], &tmp) != TCL_OK) {
	ret_val = TCL_ERROR;
	sprintf(interp -> result, "%s: non-numeric in matrix", fctn);
      } else {
	mat -> mat[j][i] = tmp;  // Matrix4 is transpose to Tcl's matrix
      }
    }
  }
  free(data_row[0]);
  free(data_row[1]);
  free(data_row[2]);
  free(data_row[3]);
  return ret_val;
}

// append the matrix into the -> result field of the interp
void tcl_append_matrix(Tcl_Interp *interp, const Matrix4 &mat)
{
  char s[TCL_DOUBLE_SPACE];

  for (int i=0; i<4; i++) {
    Tcl_AppendResult(interp, "{", NULL);
    for (int j=0; j<4; j++) {
      Tcl_PrintDouble(interp, mat.mat[j][i], s);
      Tcl_AppendResult(interp, s, (j != 3 ? " " : ""), NULL);
    }
    Tcl_AppendResult(interp, (i != 3 ? "} " : "}"), NULL);
  }
}

*/

// Given a string with a vector in it, get the vector
// YOU must delete [] the vector (in "result") when finished
// returns TCL_OK if good
// If bad, returns TCL_ERROR and sets the interp->result to the error message
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
      sprintf(interp->result, "%s: non-numeric in vector", fctn);
      ret_val = TCL_ERROR;
    } else {
      (*result)[i] = tmp;
    }
  }
  free(data);
  if (ret_val == TCL_ERROR) {
    delete [] (*result);
    *result = NULL;
  }
  return ret_val;
}

/*

// speed up the matrix * vector routines -- DIFFERENT ERROR MESSAGES
// THAN THE TCL VERSION
// speedup is nearly 25 fold
int proc_vectrans(ClientData, Tcl_Interp *interp, int argc, 
		  char *argv[])
{
  if (argc == 1) {
    Tcl_AppendResult(interp, "no value given for parameter \"m\" to \"",
		     argv[0], "\"", NULL);
    return TCL_ERROR;
  }
  if (argc == 2) {
    Tcl_AppendResult(interp, "no value given for parameter \"v\" to \"",
		     argv[0], "\"", NULL);
    return TCL_ERROR;
  }
  if (argc > 3) {
    Tcl_AppendResult(interp, "called \"", argv[0], 
		     "\" with too many arguments", NULL);
    return TCL_ERROR;
  }

  // get the matrix data
  Matrix4 mat;
  if (tcl_get_matrix(argv[0], interp, argv[1], &mat) != TCL_OK) {
    return TCL_ERROR;
  }
  // for the vector
  float *vec;
  int vec_size;
  if (tcl_get_vector(argv[0], interp, argv[2], &vec_size, 
		     &vec) != TCL_OK) {
    return TCL_ERROR;
  }
  float vec_data[4];
  if (vec_size == 3) {
    memcpy(vec_data, vec, 3*sizeof(float));
    if (!strcmp(argv[0], "coordtrans")) {
      vec_data[3] = 1;
    } else {
      vec_data[3] = 0;
    }
  } else {
    if (vec_size == 4) {
      memcpy(vec_data, vec, 4*sizeof(float));
    } else {
      Tcl_AppendResult(interp, argv[0], ": vector must be of size 3 or 4");
      delete [] vec;
      return TCL_ERROR;
    }
  }
  delete [] vec;

  // vector data is in vec_data
  float result[4];
  mat.multpoint4d(vec_data, result);
  // return it
  if (vec_size == 3) {
    char s[TCL_DOUBLE_SPACE];
    Tcl_PrintDouble(interp, result[0], s);
    Tcl_AppendElement(interp, s);
    Tcl_PrintDouble(interp, result[1], s);
    Tcl_AppendElement(interp, s);
    Tcl_PrintDouble(interp, result[2], s);
    Tcl_AppendElement(interp, s);
  } else {
    char s[TCL_DOUBLE_SPACE];
    Tcl_PrintDouble(interp, result[0], s);
    Tcl_AppendElement(interp, s);
    Tcl_PrintDouble(interp, result[1], s);
    Tcl_AppendElement(interp, s);
    Tcl_PrintDouble(interp, result[2], s);
    Tcl_AppendElement(interp, s);
    Tcl_PrintDouble(interp, result[3], s);
    Tcl_AppendElement(interp, s);
  }
  return TCL_OK;
}

// Function: transmult m1 m2 ... mn
//  Returns: the product of the matricies
// speedup is 136347 / 1316 = factor of 104
int proc_transmult(ClientData, Tcl_Interp *interp, int argc, 
		   char *argv[])
{
  // make there there are at least two values
  if (argc <= 1) {
    interp -> result = "no value given for parameter \"mx\" to \"transmult\"";
    return TCL_ERROR;
  }
  if (argc == 2) {
    interp -> result = "no value given for parameter \"my\" to \"transmult\"";
    return TCL_ERROR;
  }
  // Get the first matrix
  Matrix4 mult;
  if (tcl_get_matrix("transmult: ", interp, argv[1], &mult) != TCL_OK) {
    return TCL_ERROR;
  }
  int i = 2;
  Matrix4 tmp;
  while (i < argc) {
    if (tcl_get_matrix("transmult: ", interp, argv[i], &tmp) != TCL_OK) {
      return TCL_ERROR;
    }
    mult.multmatrix(tmp);
    i++;
  }
  tcl_append_matrix(interp, mult);
  return TCL_OK;
}

*/

#endif
