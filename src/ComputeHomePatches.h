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

#ifndef COMPUTEHOMEPATCHES_H
#define COMPUTEHOMEPATCHES_H

#include "NamdTypes.h"
#include "common.h"
#include "Compute.h"
#include "Patch.h"

#include "Box.h"
#include "OwnerBox.h"
#include "PositionBox.h"
#include "PositionOwnerBox.h"

class PatchElem {
  public:
    PatchID patchID;
    Patch *p;
    PositionBox<Patch> *positionBox;
    Box<Patch,Results> *forceBox;
    Box<Patch,AtomProperties> *atomBox;
    Position *x;
    Results *r;
    Force *f;
    AtomProperties *a;

  PatchElem() {
    patchID = -1;
    p = NULL;
    positionBox = NULL;
    forceBox = NULL;
    atomBox = NULL;
    x = NULL;
    r = NULL;
    f = NULL;
    a = NULL;
  }

  PatchElem(PatchID p) {
    patchID = p;
  }

  PatchElem(Patch *p, ComputeID cid) {
    patchID = p->getPatchID();
    this->p = p;
    positionBox = p->registerPositionPickup(cid);
    forceBox = p->registerForceDeposit(cid);
    atomBox = p->registerAtomPickup(cid);
    x = NULL;
    r = NULL;
    f = NULL;
    a = NULL;
  }
    
  ~PatchElem() {};

  int operator==(const PatchElem &a) const {
    return (a.patchID == patchID);
  }

  int operator<(const PatchElem &a) const {
    return (patchID < a.patchID);
  }
};

typedef UniqueSortedArray<PatchElem> ComputeHomePatchList;

class ReductionMgr;

class ComputeHomePatches : public Compute {
protected:
  ComputeHomePatchList patchList;

  PatchMap *patchMap;
  ReductionMgr *reduction;

public:
  ComputeHomePatches(ComputeID c);
  virtual ~ComputeHomePatches();
  virtual void initialize();
  virtual void atomUpdate();
};

#endif
/***************************************************************************
 * RCS INFORMATION:
 *
 *	$RCSfile: ComputeHomePatches.h,v $
 *	$Author: milind $	$Locker:  $		$State: Exp $
 *	$Revision: 1.1003 $	$Date: 1997/12/26 23:10:44 $
 *
 ***************************************************************************
 * REVISION HISTORY:
 *
 * $Log: ComputeHomePatches.h,v $
 * Revision 1.1003  1997/12/26 23:10:44  milind
 * Made namd2 to compile, link and run under linux. Merged Templates and src
 * directoriies, and removed separate definition and declaration files for
 * templates.
 *
 * Revision 1.1002  1997/03/18 21:35:23  jim
 * Eliminated fake_seq.  Reductions now use Patch::flags.seq.
 *
 * Revision 1.1001  1997/03/12 22:06:35  jim
 * First step towards multiple force returns and multiple time stepping.
 *
 * Revision 1.1000  1997/02/06 15:57:58  ari
 * Resetting CVS to merge branches back into the main trunk.
 * We will stick to main trunk development as suggested by CVS manual.
 * We will set up tags to track fixed points of development/release
 * as suggested by CVS manual - all praise the CVS manual.
 *
 * Revision 1.779  1997/02/06 15:52:56  ari
 * Updating Revision Line, getting rid of branches
 *
 * Revision 1.778.2.1  1997/02/05 22:18:02  ari
 * Added migration code - Currently the framework is
 * there with compiling code.  This version does
 * crash shortly after migration is complete.
 * Migration appears to complete, but Patches do
 * not appear to be left in a correct state.
 *
 * Revision 1.778  1997/01/28 01:04:20  ari
 * uplevel
 *
 * Revision 1.2  1997/01/28 00:51:47  ari
 * uplevel
 *
 * Revision 1.1.2.1  1997/01/27 21:11:40  jim
 * test
 *
 * Revision 1.777  1997/01/17 19:35:45  ari
 * Internal CVS leveling release.  Start development code work
 * at 1.777.1.1.
 *
 * Revision 1.10  1997/01/16 00:55:56  jim
 * Added reduction of energies from ComputeHomeTuples objects, except
 * for ComputeNonbondedExcl which only reports 0 energy.
 * Some problems with ReductionMgr are apparent, but it still runs.
 *
 * Revision 1.9  1997/01/14 17:59:48  jim
 * fixed multiple force additions (no adding to proxies now)
 *
 * Revision 1.8  1996/12/04 18:03:12  jim
 * added AtomProperties checkout
 *
 * Revision 1.7  1996/11/19 06:58:37  jim
 * first compiling templated version, needed ugly void* hack
 *
 * Revision 1.6  1996/11/19 04:24:24  jim
 * first templated version as ComputeHomeTuples<T>
 *
 * Revision 1.5  1996/11/18 21:28:48  ari
 * *** empty log message ***
 *
 * Revision 1.4  1996/11/04 20:06:17  nealk
 * Now it compiles :-)
 *
 * Revision 1.3  1996/11/04 19:29:02  nealk
 * Added angleForce() to system, but it is untested.
 *
 * Revision 1.2  1996/11/04 16:55:46  ari
 * *** empty log message ***
 *
 * Revision 1.1  1996/11/01 21:20:45  ari
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

