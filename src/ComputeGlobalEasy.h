/**
***  Copyright (c) 1995, 1996, 1997, 1998, 1999, 2000 by
***  The Board of Trustees of the University of Illinois.
***  All rights reserved.
**/

/*
   Forwards atoms to master node for force evaluation.
*/

#ifndef COMPUTEGLOBALEASY_H
#define COMPUTEGLOBALEASY_H

#include "ComputeGlobalMaster.h"
#include "NamdTypes.h"

class ComputeGlobalConfigMsg;
class ComputeGlobalDataMsg;
class ComputeGlobalResultsMsg;
class ComputeGlobalMaster;
class ComputeMgr;
class Molecule;
class SimParameters;
class SubmitReduction;

class ComputeGlobalEasy : public ComputeGlobalMaster {
protected:
  friend class ComputeGlobal;
  ComputeGlobalEasy(ComputeGlobal *, const char *);
  virtual ~ComputeGlobalEasy();

  int getAtomID(const char *segid, int resid, const char *aname);
  int getNumAtoms(const char* segid, int resid); // 0 on error
  int getAtomID(const char *segid, int resid, int index);
  double getMass(int atomid);
  int requestAtom(int atomid);
  int getPosition(int atomid, Position &position);
  int addForce(int atomid, Force force);
  void addEnergy(BigReal);

  virtual void easy_init(const char *);
  virtual void easy_calc(void);

private:

  virtual void initialize();
  virtual void calculate();

  ComputeGlobalConfigMsg *configMsg;
  ComputeGlobalResultsMsg *resultsMsg;
  Molecule *molecule;
  SimParameters *simParams;
  SubmitReduction *reduction;

  char *configName;
  BigReal energy;

};

#endif

