/**
***  Copyright (c) 1995, 1996, 1997, 1998, 1999, 2000 by
***  The Board of Trustees of the University of Illinois.
***  All rights reserved.
**/

/*
   HBondPairData simply stores data about a single pair of hydrogen bond
   parameters.  The atom type names can contain wildcard characters.
*/

#ifndef HBONDPAIRDATA_H
#define HBONDPAIRDATA_H

#include "common.h"

// forward references
class HBondPairData;


////////////////////////// definition of HBondPairData
class HBondPairData {

private:
  // types of atoms involved in the pair
  char *atomName[2];

  // Emin, Rmin values for the hydrogen bond
  Real Emin, Rmin;

  // computed A, B values used in calculation of energy ... these are
  // not filled in until the simulation gets running
  Real Aval, Bval;

  // flag indicating whether A and B have been yet computed
  Bool calculatedAB;

public:
  // next item in linked list, if necessary
  HBondPairData *next;

  // constructor
  HBondPairData(char *n1, char *n2, Real e, Real r) {
    atomName[0] = NAMD_stringdup(n1);
    atomName[1] = NAMD_stringdup(n2);
    Emin = e;
    Rmin = r;
    Aval = Bval = 0.0;
    calculatedAB = FALSE;
    next = NULL;
  }

  // destructor
  ~HBondPairData(void) {
    if(atomName[0])  delete [] atomName[0];
    if(atomName[1])  delete [] atomName[1];
  }

  // return data
  char *name(int i) { return atomName[i]; }
  Real emin(void) { return Emin; }
  Real rmin(void) { return Rmin; }
  Real A(void) { return Aval; }
  Real B(void) { return Bval; }
  Bool calculated(void) { return calculatedAB; }

  // save new A, B, Emin, or Rmin values
  void set_emin(Real a) { Emin = a; calculatedAB = FALSE; }
  void set_rmin(Real a) { Rmin = a; calculatedAB = FALSE; }
  void set_A(Real a) { Aval = a; calculatedAB = TRUE; }
  void set_B(Real b) { Bval = b; calculatedAB = TRUE; }
};


#endif

