//-*-c++-*-
/***************************************************************************/
/*         (C) Copyright 1996,1997 The Board of Trustees of the            */
/*                          University of Illinois                         */
/*                           All Rights Reserved                           */
/***************************************************************************/
/***************************************************************************
 * DESCRIPTION: Base class for all Compute objects. Almost an abstract
 *              class except that it does do the basic patchReady()
 *		countdown.
 *
 ***************************************************************************/

#ifndef COMPUTE_H
#define COMPUTE_H

#include "main.h"
#include "charm++.h"

#include "NamdTypes.h"

class Node;
class PatchMap;

// Base class for various forms of Compute objects
// for example: <linkto class=ComputeAngles>ComputeAngles</linkto> 
// and <linkto class=ComputeNonbondedExcl>ComputeNonbondedExcl</linkto>
class Compute {
private:
  int patchReadyCounter;
  int numPatches;
  int doAtomUpdate;
  int computeType;

protected:
  static Node* node;
  static PatchMap *patchMap;
  unsigned int basePriority;
  void enqueueWork();

public:
  const ComputeID cid;
  static int totalComputes;
  Compute(ComputeID);
  int type() { return computeType; };

  virtual ~Compute() { totalComputes--; }

  static void setNode(Node *n) { node = n; }

  void setNumPatches(int n) { patchReadyCounter = numPatches = n; }
  int getNumPatches() { return (numPatches); };

  // registers for boxes
  virtual void initialize() {};
  // destructor better unregister for boxes!

  virtual void atomUpdate() {};
  // virtual void patchReady(void);
  virtual void patchReady(PatchID pid) { if (pid > -1) patchReady(pid,0); }
  virtual void patchReady(PatchID, int);
  virtual int noWork(); // cleans up and returns 1 if no work to do
  virtual void doWork(); // actually does the work if noWork() returns 0
  virtual int sequence(void); // returns sequence number for analysis
  virtual unsigned int priority(void) { return basePriority; }
};

#endif

