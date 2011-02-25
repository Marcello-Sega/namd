/**
***  Copyright (c) 1995, 1996, 1997, 1998, 1999, 2000 by
***  The Board of Trustees of the University of Illinois.
***  All rights reserved.
**/

#ifndef GLOBALMASTERSYMMETRY_H
#define GLOBALMASTERSYMMETRY_H

#include "GlobalMaster.h"
#include "PDBData.h"
#include <map>
#include <vector>
#include "Matrix4Symmetry.h"

class GlobalMasterSymmetry : public GlobalMaster {
public:
  GlobalMasterSymmetry();
  ~GlobalMasterSymmetry();

private:
//  map <int, vector <Matrix4Symmetry> > matrices;
vector <Matrix4Symmetry> matrices;
  map < int, Matrix4Symmetry > backmatrices;

  bool gluInvertMatrix(const BigReal [16], BigReal invOut[16]);
  void initialTransform();
  void backTransform();
  void alignMonomers();
  void determineAverage();  
  void parseMatrix(int, char fileName []);
  void calculate();
  void parseAtoms(const char *file, int);

  map <int, vector < BigReal * > > averagePos;
  map <int, BigReal *> backavg;
  map <int, vector<int>  > dmap;
  map <int, BigReal * > posmap;
  map <int, BigReal * > startmap;
  map <int, BigReal > kmap;
  map <int, vector <int> > simmap;
  map <int, int> bmap;


  BigReal K;
  const char *symmetrykfile;
  int currentStep, firstStep, lastStep, firstFullStep, lastFullStep;
  bool scaleForces;
};
#endif
