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

#ifndef COMPUTE_H
#define COMPUTE_H

#include "main.h"
#include "ckdefs.h"
#include "chare.h"
#include "c++interface.h"

#include "NamdTypes.h"

#include "Templates/Box.h"
#include "Templates/OwnerBox.h"

class Patch;
class Node;
class PatchMap;

class PatchDeposit {
public:
   PatchID pid;
   Box<Patch,Force> *box;
   Patch *p;
   Force *f;

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
   Box<Patch,Position> *box;
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

// Base class for various forms of Compute objects
// including: <linkto class=ComputeAtoms>ComputeAtoms</linkto> 
// and <linkto class=ComputePatches>ComputePatches</linkto>
class Compute {
private:
  static Node* node;
  ComputeID cid;

  int numPatches;
  int patchReadyCounter;

  PatchDepositList patchDepositList;
  PatchPickupList patchPickupList;

protected :

  void enqueueWork();
  void depositAllForces();
  void registerForceDeposit(PatchID pid);
  void unregisterForceDeposit(PatchID pid);
  void registerPositionPickup(PatchID pid);
  void unregisterPositionPickup(PatchID pid);

public:
  Compute(ComputeID c) 
      : cid(c) {
  };

  ~Compute();

  static void setNode(Node *n) { node = n; }

  virtual void patchReady(void);
  virtual void patchReady(PatchID pid) { if (pid > -1) patchReady(); }
  virtual void doWork();
};

#endif
/***************************************************************************
 * RCS INFORMATION:
 *
 *	$RCSfile: Compute.h,v $
 *	$Author: ari $	$Locker:  $		$State: Exp $
 *	$Revision: 1.3 $	$Date: 1996/10/16 08:22:39 $
 *
 ***************************************************************************
 * REVISION HISTORY:
 *
 * $Log: Compute.h,v $
 * Revision 1.3  1996/10/16 08:22:39  ari
 * *** empty log message ***
 *
 * Revision 1.1  1996/08/19 22:07:49  ari
 * Initial revision
 *
 * Revision 1.4  1996/07/16 01:54:12  ari
 * *** empty log message ***
 *
 * Revision 1.3  96/07/16  01:10:26  01:10:26  ari (Aritomo Shinozaki)
 * Fixed comments, added methods
 * 
 * Revision 1.2  1996/06/25 21:10:48  gursoy
 * *** empty log message ***
 *
 * Revision 1.1  1996/06/24 14:12:26  gursoy
 * Initial revision
 *
 ***************************************************************************/

