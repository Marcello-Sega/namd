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
#include "ckdefs.h"
#include "chare.h"
#include "c++interface.h"

#include "NamdTypes.h"
#include "Compute.h"

#include "Templates/Box.h"
#include "Templates/OwnerBox.h"
#include "PositionBox.h"
#include "PositionOwnerBox.h"

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
/***************************************************************************
 * RCS INFORMATION:
 *
 *	$RCSfile: ComputeGeneral.h,v $
 *	$Author: ari $	$Locker:  $		$State: Exp $
 *	$Revision: 1.778 $	$Date: 1997/01/28 00:30:09 $
 *
 ***************************************************************************
 * REVISION HISTORY:
 *
 * $Log: ComputeGeneral.h,v $
 * Revision 1.778  1997/01/28 00:30:09  ari
 * internal release uplevel to 1.778
 *
 * Revision 1.777.2.2  1997/01/27 22:44:59  ari
 * Basic Atom Migration Code added.
 * Added correct magic first line to .h files for xemacs to go to C++ mode.
 * Compiles and runs without migration turned on.
 *
 * Revision 1.777.2.1  1997/01/24 22:00:26  jim
 * Changes for periodic boundary conditions.
 *
 * Revision 1.777  1997/01/17 19:35:42  ari
 * Internal CVS leveling release.  Start development code work
 * at 1.777.1.1.
 *
 * Revision 1.1  1996/10/22 19:15:14  ari
 * Initial revision
 *
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

