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

#ifndef COMPUTEHOMETUPLES_H
#define COMPUTEHOMETUPLES_H

#include "NamdTypes.h"
#include "common.h"
#include "structures.h"
#include "Compute.h"
#include "Patch.h"

#include "Templates/Box.h"
#include "Templates/OwnerBox.h"
#include "PositionBox.h"
#include "PositionOwnerBox.h"

enum PatchType {HOME,PROXY};

class TuplePatchElem {
  public:
    PatchID patchID;
    Patch *p;
    PatchType patchType;
    PositionBox<Patch> *positionBox;
    Box<Patch,Force> *forceBox;
    Box<Patch,AtomProperties> *atomBox;
    Position *x;
    Force *f;
    AtomProperties *a;

  TuplePatchElem() {
    patchID = -1;
    p = NULL;
    positionBox = NULL;
    forceBox = NULL;
    atomBox = NULL;
    x = NULL;
    f = NULL;
    a = NULL;
  }

  TuplePatchElem(PatchID p) {
    patchID = p;
  }

  TuplePatchElem(Patch *p, PatchType pt, ComputeID cid) {
    patchID = p->getPatchID();
    this->p = p;
    patchType = pt;
    positionBox = p->registerPositionPickup(cid);
    if ( pt == HOME ) forceBox = p->registerForceDeposit(cid);
    else forceBox =  NULL;
    atomBox = p->registerAtomPickup(cid);
    x = NULL;
    f = NULL;
    a = NULL;
  }
    
  ~TuplePatchElem() {};

  int operator==(const TuplePatchElem &a) const {
    return (a.patchID == patchID);
  }

  int operator<(const TuplePatchElem &a) const {
    return (patchID < a.patchID);
  }
};

typedef UniqueSortedArray<TuplePatchElem> TuplePatchList;

class AtomMap;
class ReductionMgr;

template <class T>
class ComputeHomeTuples : public Compute {
private:
  void loadTuples();
  void sizeDummy();

  UniqueSortedArray<T> tupleList;
  TuplePatchList tuplePatchList;

  PatchMap *patchMap;
  AtomMap *atomMap;
  ReductionMgr *reduction;

  int fake_seq;

  int maxProxyAtoms;
  Force *dummy;
  
public:
  ComputeHomeTuples(ComputeID c);
  virtual ~ComputeHomeTuples();
  void initialize();
  void atomUpdate();
  void doWork();
};

#endif
/***************************************************************************
 * RCS INFORMATION:
 *
 *	$RCSfile: ComputeHomeTuples.h,v $
 *	$Author: ari $	$Locker:  $		$State: Exp $
 *	$Revision: 1.1000 $	$Date: 1997/02/06 15:58:00 $
 *
 ***************************************************************************
 * REVISION HISTORY:
 *
 * $Log: ComputeHomeTuples.h,v $
 * Revision 1.1000  1997/02/06 15:58:00  ari
 * Resetting CVS to merge branches back into the main trunk.
 * We will stick to main trunk development as suggested by CVS manual.
 * We will set up tags to track fixed points of development/release
 * as suggested by CVS manual - all praise the CVS manual.
 *
 * Revision 1.779  1997/02/06 15:52:58  ari
 * Updating Revision Line, getting rid of branches
 *
 * Revision 1.778.2.1  1997/02/05 22:18:04  ari
 * Added migration code - Currently the framework is
 * there with compiling code.  This version does
 * crash shortly after migration is complete.
 * Migration appears to complete, but Patches do
 * not appear to be left in a correct state.
 *
 * Revision 1.778  1997/01/28 00:30:12  ari
 * internal release uplevel to 1.778
 *
 * Revision 1.777.2.2  1997/01/27 22:45:01  ari
 * Basic Atom Migration Code added.
 * Added correct magic first line to .h files for xemacs to go to C++ mode.
 * Compiles and runs without migration turned on.
 *
 * Revision 1.777.2.1  1997/01/24 22:00:28  jim
 * Changes for periodic boundary conditions.
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

