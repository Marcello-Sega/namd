//-*-c++-*-
/***************************************************************************/
/*                                                                         */
/*              (C) Copyright 1996 The Board of Trustees of the            */
/*                          University of Illinois                         */
/*                           All Rights Reserved                           */
/*									   */
/***************************************************************************/

/***************************************************************************
 * DESCRIPTION:
 *
 ***************************************************************************/

#ifndef COMPUTEGENERAL_H
#define COMPUTEGENERAL_H

#include "main.h"
#include "charm++.h"

#include "NamdTypes.h"
#include "PatchTypes.h"
#include "Compute.h"

#include "Box.h"
#include "OwnerBox.h"
#include "PositionBox.h"
#include "PositionOwnerBox.h"

class Patch;
class Node;
class PatchMap;

class PatchDeposit {
public:
   PatchID pid;
   Box<Patch,Results> *box;
   Patch *p;
   Results *r;

   PatchDeposit(PatchID p) : pid(p) {};
   PatchDeposit() : pid(-1) {};
   ~PatchDeposit() {};
   int operator== (const PatchDeposit p) const { 
     return (pid == p.pid);
   }
   int operator< (const PatchDeposit p) const {
     return (pid < p.pid);
   }
};

class PatchPickup {
public:
   PatchID pid;
   PositionBox<Patch> *box;
   Patch *p;
   Position *x;

   PatchPickup(PatchID p) : pid(p) {};
   PatchPickup() : pid(-1) {};
   ~PatchPickup() {};

   int operator== (const PatchPickup p) const { 
     return (pid == p.pid);
   }
   int operator< (const PatchPickup p) const {
     return (pid < p.pid);
   }
};

typedef ResizeArrayIter<PatchDeposit> PatchDepositListIter;
typedef SortedArray<PatchDeposit> PatchDepositList;
typedef ResizeArrayIter<PatchPickup> PatchPickupListIter;
typedef SortedArray<PatchPickup> PatchPickupList;

class ComputeGeneral : public Compute {
private:
  PatchDepositList patchDepositList;
  PatchPickupList patchPickupList;

protected :
  void depositAllForces();
  void registerForceDeposit(PatchID pid);
  void unregisterForceDeposit(PatchID pid);
  void registerPositionPickup(PatchID pid);
  void unregisterPositionPickup(PatchID pid);

public:
  ComputeGeneral(ComputeID c) : Compute(c) {};
  virtual ~ComputeGeneral();

  virtual void doWork();
};

#endif

