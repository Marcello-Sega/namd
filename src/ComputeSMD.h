/**
***  Copyright (c) 1995, 1996, 1997, 1998, 1999, 2000 by
***  The Board of Trustees of the University of Illinois.
***  All rights reserved.
**/

#ifndef COMPUTESMD_H
#define COMPUTESMD_H

#include "ComputeGlobalMaster.h"

class ComputeGlobalConfigMsg;
class ComputeGlobalResultsMsg;

class ComputeSMD : public ComputeGlobalMaster {
public:
  ComputeSMD(ComputeMgr *);
  ~ComputeSMD();

private:
  virtual void initialize();
  virtual void calculate();

  void output(int, Position, Force);
  void parse_atoms(char *);
 
  ComputeGlobalConfigMsg *configMsg;

  BigReal k;
  BigReal moveVel;   // A/timestep
  Vector moveDir;
  int outputFreq;
  Position cm;       // Initial center of mass
  int currentTime;   // Keep track of elapsed time steps for yourself!
};
#endif

