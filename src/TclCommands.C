/**
***  Copyright (c) 1995, 1996, 1997, 1998, 1999, 2000 by
***  The Board of Trustees of the University of Illinois.
***  All rights reserved.
**/

#include <stdlib.h> 
#ifndef _NO_MALLOC_H
#include <malloc.h>
#endif
#include <errno.h>
#include "TclCommands.h"

#ifdef NAMD_TCL

#include <tcl.h>
#include <Vector.h>
#include <NamdTypes.h>

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


// Get a 3-D vector from a TCL list
int get_3D_vector(Tcl_Interp *interp, char *list, Vector &result)
{
  int num, status=1;
  char **data;
  
  if (Tcl_SplitList(interp,list,&num,&data) != TCL_OK)
    status = 0;
  else if (num != 3)
    status = 0;
  else if (Tcl_GetDouble(interp,data[0],&(result.x)) != TCL_OK)
    status = 0;
  else if (Tcl_GetDouble(interp,data[1],&(result.y)) != TCL_OK)
    status = 0;
  else if (Tcl_GetDouble(interp,data[2],&(result.z)) != TCL_OK)
    status = 0;
  
  Tcl_Free((char*) data);
  return status;
}


// Append a 3-D vector to the result string
void append_3D_vector(Tcl_Interp *interp, Vector &v)
{
  char s[3][TCL_DOUBLE_SPACE], *t[3], *list;
  
  Tcl_PrintDouble(interp,v.x,s[0]);
  Tcl_PrintDouble(interp,v.y,s[1]);
  Tcl_PrintDouble(interp,v.z,s[2]);
  t[0]=s[0], t[1]=s[1], t[2]=s[2];
  list = Tcl_Merge(3,t);
  Tcl_AppendElement(interp,list);
  Tcl_Free(list);
}


// Function: getbond coor1 coor2
//  Returns: the length of the bond formed by the two atoms (i.e., the distance between them)
int proc_getbond(ClientData, Tcl_Interp *interp, int argc, char *argv[])
{
  char s[TCL_DOUBLE_SPACE];
  Vector r1, r2;
  
  if (argc != 3)
    return TCL_ERROR;
  if (!get_3D_vector(interp,argv[1],r1))
    return TCL_ERROR;
  if (!get_3D_vector(interp,argv[2],r2))
    return TCL_ERROR;
  Tcl_PrintDouble(interp,(r2-r1).length(),s);
  Tcl_SetResult(interp,s,TCL_VOLATILE);
  return TCL_OK;
}


// Function: getangle coor1 coor2 coor3
//  Returns: the angle formed by the three atoms
int proc_getangle(ClientData, Tcl_Interp *interp, int argc, char *argv[])
{
  char s[TCL_DOUBLE_SPACE];
  Vector r1, r2, r3, r12, r32;
  
  if (argc != 4)
    return TCL_ERROR;
  if (!get_3D_vector(interp,argv[1],r1))
    return TCL_ERROR;
  if (!get_3D_vector(interp,argv[2],r2))
    return TCL_ERROR;
  if (!get_3D_vector(interp,argv[3],r3))
    return TCL_ERROR;
  r12 = r1 - r2;
  r32 = r3 - r2;
  Tcl_PrintDouble(interp,acos((r12*r32)/(r12.length()*r32.length()))*180/PI,s);
  Tcl_SetResult(interp,s,TCL_VOLATILE);
  return TCL_OK;
}


// Function: getdihedral coor1 coor2 coor3 coor4
//  Returns: the dihedral formed by the four atoms
int proc_getdihedral(ClientData, Tcl_Interp *interp, int argc, char *argv[])
{
  char s[TCL_DOUBLE_SPACE];
  BigReal rA, rB, rC;
  Vector r1, r2, r3, r4, r12, r23, r34, A, B, C;
  
  if (argc != 5)
    return TCL_ERROR;
  if (!get_3D_vector(interp,argv[1],r1))
    return TCL_ERROR;
  if (!get_3D_vector(interp,argv[2],r2))
    return TCL_ERROR;
  if (!get_3D_vector(interp,argv[3],r3))
    return TCL_ERROR;
  if (!get_3D_vector(interp,argv[4],r4))
    return TCL_ERROR;
  r12 = r1 - r2;
  r23 = r2 - r3;
  r34 = r3 - r4;
  A = cross(r12,r23);
  B = cross(r23,r34);
  C = cross(r23,A);
  rA = A.length();
  rB = B.length();
  rC = C.length();
  Tcl_PrintDouble(interp,-atan2((C*B)/(rC*rB),(A*B)/(rA*rB))*180/PI,s);
  Tcl_SetResult(interp,s,TCL_VOLATILE);
  return TCL_OK;
}


// Function: anglegrad coor1 coor2 coor3
//  Returns: a list of gradients for each atom
// The code was basically copied from ComputeAngles.C
int proc_anglegrad(ClientData, Tcl_Interp *interp, int argc, char *argv[])
{
  Vector r1, r2, r3, r12, r32;
  
  if (argc != 4)
    return TCL_ERROR;
  if (!get_3D_vector(interp,argv[1],r1))
    return TCL_ERROR;
  if (!get_3D_vector(interp,argv[2],r2))
    return TCL_ERROR;
  if (!get_3D_vector(interp,argv[3],r3))
    return TCL_ERROR;
  
  r12 = r1 - r2;
  BigReal d12 = r12.length();
  r32 = r3 - r2;
  BigReal d32 = r32.length();
  
  BigReal cos_theta = (r12*r32)/(d12*d32);
  
  //  Normalize vector r12 and r32
  BigReal d12inv = 1. / d12;
  BigReal d32inv = 1. / d32;
  
  BigReal sin_theta = sqrt(1.0 - cos_theta*cos_theta);
  BigReal diff = 1/sin_theta;
  BigReal c1 = diff * d12inv;
  BigReal c2 = diff * d32inv;
  
  //  Calculate the actual forces
  Force force1 = c1*(r12*(d12inv*cos_theta) - r32*d32inv);
  Force force2 = force1;
  Force force3 = c2*(r32*(d32inv*cos_theta) - r12*d12inv);
  force2 += force3;  force2 *= -1;
  
  append_3D_vector(interp,force1);
  append_3D_vector(interp,force2);
  append_3D_vector(interp,force3);
  
  return TCL_OK;
}


// Function: dihedralgrad coor1 coor2 coor3 coor4
//  Returns: a list of gradients for each atom
// The code was basically copied from ComputeDihedrals.C
int proc_dihedralgrad(ClientData, Tcl_Interp *interp, int argc, char *argv[])
{
  BigReal K1;
  Vector r1, r2, r3, r4, r12, r23, r34;
  
  if (argc != 5)
    return TCL_ERROR;
  if (!get_3D_vector(interp,argv[1],r1))
    return TCL_ERROR;
  if (!get_3D_vector(interp,argv[2],r2))
    return TCL_ERROR;
  if (!get_3D_vector(interp,argv[3],r3))
    return TCL_ERROR;
  if (!get_3D_vector(interp,argv[4],r4))
    return TCL_ERROR;
  
  r12 = r1 - r2;
  r23 = r2 - r3;
  r34 = r3 - r4;
  
  //  Calculate the cross products and distances
  Vector A = cross(r12,r23);
  BigReal rA = A.length();
  Vector B = cross(r23,r34);
  BigReal rB = B.length();
  Vector C = cross(r23,A);
  BigReal rC = C.length();
  
  //  Calculate the sin and cos
  BigReal cos_phi = (A*B)/(rA*rB);
  BigReal sin_phi = (C*B)/(rC*rB);
  
  Force f1,f2,f3;
  
  //  Normalize B
  rB = 1.0/rB;
  B *= rB;
  
  //  We first need to figure out whether the
  //  sin or cos form will be more stable.  For this,
  //  just look at the value of phi
  if (fabs(sin_phi) > 0.1)
  {
    //  use the sin version to avoid 1/cos terms
    
    //  Normalize A
    rA = 1.0/rA;
    A *= rA;
    Vector dcosdA = rA*(cos_phi*A-B);
    Vector dcosdB = rB*(cos_phi*B-A);
    
    K1 = -1/sin_phi;
    
    f1.x = K1*(r23.y*dcosdA.z - r23.z*dcosdA.y);
    f1.y = K1*(r23.z*dcosdA.x - r23.x*dcosdA.z);
    f1.z = K1*(r23.x*dcosdA.y - r23.y*dcosdA.x);
    
    f3.x = K1*(r23.z*dcosdB.y - r23.y*dcosdB.z);
    f3.y = K1*(r23.x*dcosdB.z - r23.z*dcosdB.x);
    f3.z = K1*(r23.y*dcosdB.x - r23.x*dcosdB.y);
    
    f2.x = K1*(r12.z*dcosdA.y - r12.y*dcosdA.z
             + r34.y*dcosdB.z - r34.z*dcosdB.y);
    f2.y = K1*(r12.x*dcosdA.z - r12.z*dcosdA.x
             + r34.z*dcosdB.x - r34.x*dcosdB.z);
    f2.z = K1*(r12.y*dcosdA.x - r12.x*dcosdA.y
             + r34.x*dcosdB.y - r34.y*dcosdB.x);
  }
  else
  {
    //  This angle is closer to 0 or 180 than it is to
    //  90, so use the cos version to avoid 1/sin terms
    
    //  Normalize C
    rC = 1.0/rC;
    C *= rC;
    Vector dsindC = rC*(sin_phi*C-B);
    Vector dsindB = rB*(sin_phi*B-C);
    
    K1 = 1/cos_phi;
    
    f1.x = K1*((r23.y*r23.y + r23.z*r23.z)*dsindC.x
              - r23.x*r23.y*dsindC.y
              - r23.x*r23.z*dsindC.z);
    f1.y = K1*((r23.z*r23.z + r23.x*r23.x)*dsindC.y
              - r23.y*r23.z*dsindC.z
              - r23.y*r23.x*dsindC.x);
    f1.z = K1*((r23.x*r23.x + r23.y*r23.y)*dsindC.z
              - r23.z*r23.x*dsindC.x
              - r23.z*r23.y*dsindC.y);
    
    f3 = cross(K1,dsindB,r23);
    
    f2.x = K1*(-(r23.y*r12.y + r23.z*r12.z)*dsindC.x
           +(2.0*r23.x*r12.y - r12.x*r23.y)*dsindC.y
           +(2.0*r23.x*r12.z - r12.x*r23.z)*dsindC.z
           +dsindB.z*r34.y - dsindB.y*r34.z);
    f2.y = K1*(-(r23.z*r12.z + r23.x*r12.x)*dsindC.y
           +(2.0*r23.y*r12.z - r12.y*r23.z)*dsindC.z
           +(2.0*r23.y*r12.x - r12.y*r23.x)*dsindC.x
           +dsindB.x*r34.z - dsindB.z*r34.x);
    f2.z = K1*(-(r23.x*r12.x + r23.y*r12.y)*dsindC.z
           +(2.0*r23.z*r12.x - r12.z*r23.x)*dsindC.x
           +(2.0*r23.z*r12.y - r12.z*r23.y)*dsindC.y
           +dsindB.y*r34.x - dsindB.x*r34.y);
  }
  
  append_3D_vector(interp,f1);
  append_3D_vector(interp,f2-f1);
  append_3D_vector(interp,f3-f2);
  append_3D_vector(interp,-f3);
  
  return TCL_OK;
}

#endif
