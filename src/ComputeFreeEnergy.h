/**
***  Copyright (c) 1995, 1996, 1997, 1998, 1999, 2000 by
***  The Board of Trustees of the University of Illinois.
***  All rights reserved.
**/

/*
   Forwards atoms to master node for force evaluation.
*/

#if !defined(COMPUTEFREEENERGY_H)
  #define COMPUTEFREEENERGY_H

#include "ComputeGlobalMaster.h"
#include "NamdTypes.h"
#include "iostream.h"

class ComputeGlobalConfigMsg;
class ComputeGlobalDataMsg;
class ComputeGlobalResultsMsg;
class ComputeGlobalMaster;
class ComputeMgr;
class Molecule;
class SimParameters;

class ComputeFreeEnergy : public ComputeGlobalMaster {
public:
  ComputeFreeEnergy(ComputeMgr *);
  ~ComputeFreeEnergy();
private:
  virtual void initialize();
  virtual void calculate();
  virtual void user_initialize();
  virtual void user_calculate();
  void update();
  ComputeGlobalConfigMsg *configMsg;
  ComputeGlobalResultsMsg *resultsMsg;
  Molecule *molecule;
  SimParameters *simParams;
  char *config;
  ARestraintManager  m_RestraintManager;
  ALambdaManager     m_LambdaManager;
public:
  // These all return -1 on error.
  int getAtomID(const char *segid, int resid, const char *aname);
  int getNumAtoms(const char* segid, int resid); // 0 on error
  int getAtomID(const char *segid, int resid, int index);
  double getMass(int atomid);
  int requestAtom(int atomid);
  int getPosition(int atomid, Position &position);
  int addForce(int atomid, Force force);
};

#endif

