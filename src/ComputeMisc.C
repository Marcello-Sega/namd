/**
***  Copyright (c) 1995, 1996, 1997, 1998, 1999, 2000 by
***  The Board of Trustees of the University of Illinois.
***  All rights reserved.
**/

/*
   Forwards atoms to master node for force evaluation.
*/

#include "ComputeMisc.h"
#include "InfoStream.h"

ComputeMisc::ComputeMisc(ComputeMgr *c)
  : ComputeGlobalEasy(c,"miscForcesScript") {
  ;
}

ComputeMisc::~ComputeMisc() {
  ;
}

/* you may use the following in easy_init or easy_calc
  int getAtomID(const char *segid, int resid, const char *aname); // zero-based
  int getNumAtoms(const char* segid, int resid); // 0 on error
  int getAtomID(const char *segid, int resid, int index);
  double getMass(int atomid); // -1.0 on error
*/

/* you may use the following only in easy_init()
  int requestAtom(int atomid); // 0 on success, -1 on error
*/

void ComputeMisc::easy_init(const char *config) {
  iout << iINFO << "  MISC FORCES CONFIG\n";
  iout << iINFO << "**********************\n";
  iout << config;
  iout << iINFO << "**********************\n" << endi;

  requestAtom(0);
}

/* you may use the following only in easy_calc()
  int getPosition(int atomid, Position &position); // 0 on success, -1 on error
  int addForce(int atomid, Force force); // 0 on success, -1 on error
  void addEnergy(BigReal);
*/

void ComputeMisc::easy_calc() {
  Vector myp;
  BigReal k = 10.0;
  getPosition(0,myp);
  iout << iINFO << "Atom 0 is at " << myp << "\n" << endi;
  addForce(0,-k * myp);
  addEnergy(0.5 * k * myp.length2());
}


