/**
***  Copyright (c) 1995, 1996, 1997, 1998, 1999, 2000 by
***  The Board of Trustees of the University of Illinois.
***  All rights reserved.
**/

#ifndef GLOBALMASTERTMD_H
#define GLOBALMASTERTMD_H

#include "GlobalMaster.h"

class GlobalMasterTMD : public GlobalMaster {
public:
  GlobalMasterTMD();
  ~GlobalMasterTMD();

private:
  virtual void calculate();
  void parseAtoms(const char *file, int);

  int numTMDatoms;
  BigReal k;
  BigReal initialRMS, finalRMS;
  int outputFreq;
  int currentStep, firstStep, lastStep;
  BigReal *target;

  // mapping of atom id's to array positions
  int *aidmap;
};
#endif

