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

#ifndef COMPUTEANGLE_H
#define COMPUTEANGLE_H

#include "NamdTypes.h"
#include "common.h"
#include "structures.h"
#include "Compute.h"
#include "Patch.h"

#include "Templates/Box.h"
#include "Templates/OwnerBox.h"

enum PatchType {HOME,PROXY};

class TuplePatchElem {
  public:
    PatchID patchID;
    Patch *p;
    PatchType patchType;
    Box<Patch,Position> *positionBox;
    Box<Patch,Force> *forceBox;
    Position *x;
    Force *f;

  TuplePatchElem() {
    patchID = -1;
    p = NULL;
    positionBox = NULL;
    forceBox = NULL;
    x = NULL;
    f = NULL;
  }

  TuplePatchElem(PatchID p) {
    patchID = p;
  }

  TuplePatchElem(Patch *p, PatchType pt, ComputeID cid) {
    patchID = p->getPatchID();
    this->p = p;
    patchType = pt;
    positionBox = p->registerPositionPickup(cid);
    forceBox = p->registerForceDeposit(cid);
    x = NULL;
    f = NULL;
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

template <class T>
class ComputeHomeTuples : public Compute {
private:
  UniqueSortedArray<T> tupleList;
  TuplePatchList tuplePatchList;

  PatchMap *patchMap;
  AtomMap *atomMap;

  int maxProxyAtoms;
  Force *dummy;
  
public:
  ComputeHomeTuples(ComputeID c);
  virtual ~ComputeHomeTuples() {
    delete [] dummy;
  }

  void mapReady();
  void doWork();
};

#endif
/***************************************************************************
 * RCS INFORMATION:
 *
 *	$RCSfile: ComputeHomeTuples.h,v $
 *	$Author: jim $	$Locker:  $		$State: Exp $
 *	$Revision: 1.6 $	$Date: 1996/11/19 04:24:24 $
 *
 ***************************************************************************
 * REVISION HISTORY:
 *
 * $Log: ComputeHomeTuples.h,v $
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

