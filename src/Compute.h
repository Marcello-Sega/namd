/**
***  Copyright (c) 1995, 1996, 1997, 1998, 1999, 2000 by
***  The Board of Trustees of the University of Illinois.
***  All rights reserved.
**/

/*
   Base class for all Compute objects. Almost an abstract
   class except that it does do the basic patchReady()
   countdown.
*/

#ifndef COMPUTE_H
#define COMPUTE_H

#include "main.h"
#include "charm++.h"

#include "NamdTypes.h"

class Node;
class PatchMap;
class LocalWorkMsg;

// Base class for various forms of Compute objects
// for example: <linkto class=ComputeAngles>ComputeAngles</linkto> 
// and <linkto class=ComputeNonbondedExcl>ComputeNonbondedExcl</linkto>
class Compute {
private:
  int patchReadyCounter;
  int numPatches;
  int doAtomUpdate;
  int computeType;
  int sequenceNumber;

protected:
  static Node* node;
  static PatchMap *patchMap;
  unsigned int basePriority;
  void enqueueWork();

public:
  const ComputeID cid;
  LocalWorkMsg *const localWorkMsg;
  static int totalComputes;
  Compute(ComputeID);
  int type() { return computeType; };

  virtual ~Compute();

  static void setNode(Node *n) { node = n; }

  void setNumPatches(int n) { patchReadyCounter = numPatches = n; }
  int getNumPatches() { return (numPatches); };

  // registers for boxes
  virtual void initialize() {};
  // destructor better unregister for boxes!

  virtual void atomUpdate() {};
  virtual void patchReady(PatchID, int doneMigration, int seq);
  virtual int noWork(); // cleans up and returns 1 if no work to do
  virtual void doWork(); // actually does the work if noWork() returns 0
  int sequence(void) { return sequenceNumber; }
  virtual unsigned int priority(void) { return basePriority; }
};

#endif

