/**
***  Copyright (c) 1995, 1996, 1997, 1998, 1999, 2000 by
***  The Board of Trustees of the University of Illinois.
***  All rights reserved.
**/

#ifndef COMPUTETCLBC_H
#define COMPUTETCLBC_H

#ifdef NAMD_TCL
#include <tcl.h>
#endif

#include "ComputeHomePatches.h"
#include "ReductionMgr.h"
#include "Tensor.h"

class ComputeMgr;

class ComputeTclBC : public ComputeHomePatches {

public:
  ComputeTclBC(ComputeID c);
  virtual ~ComputeTclBC();
  void doWork();

private:
  ResizeArrayIter<PatchElem> ap;
  int i_atom, n_atom;
  CompAtom *atoms;
  FullAtom *fullatoms;
  Force *forces;
  BigReal energy;
  SubmitReduction *reduction;

#ifdef NAMD_TCL
  Tcl_Interp *interp;
  static int Tcl_print(ClientData, Tcl_Interp *, int, char **);
  static int Tcl_nextatom(ClientData, Tcl_Interp *, int, Tcl_Obj * const []);
  static int Tcl_getcoord(ClientData, Tcl_Interp *, int, Tcl_Obj * const []);
  static int Tcl_getmass(ClientData, Tcl_Interp *, int, Tcl_Obj * const []);
  static int Tcl_getcharge(ClientData, Tcl_Interp *, int, Tcl_Obj * const []);
  static int Tcl_getid(ClientData, Tcl_Interp *, int, Tcl_Obj * const []);
  static int Tcl_addforce(ClientData, Tcl_Interp *, int, Tcl_Obj * const []);
  static int Tcl_addenergy(ClientData, Tcl_Interp *, int, Tcl_Obj * const []);
#endif

};

#endif







