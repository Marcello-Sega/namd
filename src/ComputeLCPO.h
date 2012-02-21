/**
***  Copyright (c) 1995, 1996, 1997, 1998, 1999, 2000 by
***  The Board of Trustees of the University of Illinois.
***  All rights reserved.
**/

/*
   Solvent Accessible Surface Area calculation by LCPO
   Linear Combination of Pairwise Overlaps
*/

#ifndef COMPUTELCPO_H
#define COMPUTELCPO_H

#include "Compute.h"
#include "PatchTypes.h"
#include "Box.h"
#include "OwnerBox.h"
#include "ComputeNonbondedUtil.h"

class Patch;
class Node;
class PatchMap;

class ComputeLCPO: public Compute, private ComputeNonbondedUtil {

public:
  ComputeLCPO(ComputeID c, PatchID pid[], int t[],
		ComputeNonbondedWorkArrays* _workArrays,
		int minPartition, int maxPartition, int numPartitions, int numPatches);

  virtual ~ComputeLCPO();

  virtual void initialize();
  virtual void atomUpdate();
  virtual void doWork();
  virtual int noWork();

protected :
  int numAtoms[8];
  //int gbisPhase;
  CompAtomExt *posExt[8];
  CompAtom *pos[8];
  Results *force[8];
  int *lcpoType[8];

  int pairlistsMaxAge;
  int pairlistsAge;
  int step;

  virtual void doForce();
  Patch *patch[8];

  PatchID patchID[8];
  int trans[8];
  Box<Patch,CompAtom> *positionBox[8];
  Box<Patch,Results> *forceBox[8];
  Box<Patch,int> *lcpoTypeBox[8];

  ComputeNonbondedWorkArrays* const workArrays;
  int minPart, maxPart, numParts;

  BigReal reductionData[reductionDataSize];
  SubmitReduction *reduction;

  private:
    int octet[8];//maps 0-7 into patchID for invalid patches
    int pid8[8];
    BigReal bounds[3][2];
    int periodic[3];
    int oob[3];
    nonbonded params;
    Vector offset[8];
    int minIg[8];
    int strideIg;//stride through partitions

    //index "i" is patch; index "j" is valid atoms in patch
    Pairlists inAtomsPl;
    //a pairlist for each patch
    Pairlists pairlists[8];
    //Pairlists triplets[8];
    Real surfTen;

    static const Real lcpoParams[23][5];

    int isInBounds( Real x, Real y, Real z );
};

#endif
