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

#ifndef COMPUTEBOND_H
#define COMPUTEBOND_H

#include "ComputeHomeTuples.h"
class Molecule;

class BondElem {
public:
    // ComputeHomeTuples interface
    enum { size = 2 };
    AtomID atomID[size];
    int    localIndex[size];
    TuplePatchElem *p[size];
    BigReal computeForce(void);
    // The following is evil, but the compiler chokes otherwise. (JCP)
    static void addTuplesForAtom(void*, AtomID, Molecule*);

    // Internal data
    Index bondType;


  BondElem() {
    atomID[0] = -1;
    atomID[1] = -1;
    p[0] = NULL;
    p[1] = NULL;
  }
  BondElem(const Bond *a) {
    atomID[0] = a->atom1;
    atomID[1] = a->atom2;
    bondType = a->bond_type;
  }

  BondElem(AtomID atom0, AtomID atom1) {
    if (atom0 > atom1) {  // Swap end atoms so lowest is first!
      AtomID tmp = atom1; atom1 = atom0; atom0 = tmp; 
    }
    atomID[0] = atom0;
    atomID[1] = atom1;
  }
  ~BondElem() {};

  int operator==(const BondElem &a) const {
    return (a.atomID[0] == atomID[0] && a.atomID[1] == atomID[1]);
  }

  int operator<(const BondElem &a) const {
    return (atomID[0] < a.atomID[0] ||
            (atomID[0] == a.atomID[0] &&
            (atomID[1] < a.atomID[1]) ));
  }
};

class ComputeBonds : public ComputeHomeTuples<BondElem>
{
public:

  ComputeBonds(ComputeID c) : ComputeHomeTuples<BondElem>(c) { ; }

};

#endif
/***************************************************************************
 * RCS INFORMATION:
 *
 *	$RCSfile: ComputeBonds.h,v $
 *	$Author: nealk $	$Locker:  $		$State: Exp $
 *	$Revision: 1.2 $	$Date: 1996/12/03 15:15:40 $
 *
 ***************************************************************************
 * REVISION HISTORY:
 *
 * $Log: ComputeBonds.h,v $
 * Revision 1.2  1996/12/03 15:15:40  nealk
 * Removed tons-o-debugging.
 *
 * Revision 1.1  1996/12/03 14:53:42  nealk
 * Initial revision
 *
 *
 ***************************************************************************/

