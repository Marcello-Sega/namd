/**
***  Copyright (c) 1995, 1996, 1997, 1998, 1999, 2000 by
***  The Board of Trustees of the University of Illinois.
***  All rights reserved.
**/

/*
   Compute object which deals with a single patch.
*/

#ifndef COMPUTEPATCH_H
#define COMPUTEPATCH_H

#include "Compute.h"
#include "PatchTypes.h"

#include "Box.h"
#include "OwnerBox.h"
#include "PositionBox.h"
#include "PositionOwnerBox.h"

class Patch;
class HomePatch;
class Node;
class PatchMap;

class ComputePatch : public Compute {

public:
  ComputePatch(ComputeID c, PatchID pid);
  virtual ~ComputePatch();

  virtual void initialize();
  virtual void atomUpdate();
  virtual void doWork();
  virtual int sequence(void); // returns sequence number for analysis

protected :
  int numAtoms;
  virtual void doForce(Position* p, Results* r, AtomProperties* a);
  Patch *patch;

private:
  PatchID patchID;
  PositionBox<Patch> *positionBox;
  Box<Patch,Results> *forceBox;
  Box<Patch,AtomProperties> *atomBox;

friend class ComputeHomePatch;

};

class ComputeHomePatch : public ComputePatch {

public:
  ComputeHomePatch(ComputeID c, PatchID pid);
  virtual void initialize();
  virtual void doWork();

protected:
  HomePatch *homePatch;
  virtual void doForce(Position* p, Results* r, AtomProperties* a, Transform* t);

};

#endif

