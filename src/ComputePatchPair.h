//-*-c++-*-
/***************************************************************************/
/*          (C) Copyright 1996,1997 The Board of Trustees of the           */
/*                          University of Illinois                         */
/*                           All Rights Reserved                           */
/***************************************************************************/
/***************************************************************************
 * DESCRIPTION: Primary class for pairwise force calculations over
 *              pairs of patches.  Takes care of boxes, depositing of
 *              Forces etc.
 *
 ***************************************************************************/

#ifndef COMPUTEPPAIR_H
#define COMPUTEPPAIR_H

#include "Compute.h"
#include "PatchTypes.h"

#include "Box.h"
#include "OwnerBox.h"
#include "PositionBox.h"
#include "PositionOwnerBox.h"

class Patch;
class Node;
class PatchMap;

class ComputePatchPair : public Compute {

public:
  ComputePatchPair(ComputeID c, PatchID pid[], int t[]);
  virtual ~ComputePatchPair();

  virtual void initialize();
  virtual void atomUpdate();
  virtual void doWork();
  virtual int sequence(void); // returns sequence number for analysis

protected :
  int numAtoms[2];
  virtual void doForce(Position* p[2], Results* r[2], AtomProperties* a[2]);
  Patch *patch[2];

// private: // hack for ComputeNonbondedPair::noWork()
  PatchID patchID[2];
  int trans[2];
  PositionBox<Patch> *positionBox[2];
  Box<Patch,AtomProperties> *atomBox[2];
  Box<Patch,Results> *forceBox[2];
};

#endif

