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

#ifndef COMPUTEBOND_INL
#define COMPUTEBOND_INL

#include "ComputeBonds.h"

inline BondElem::BondElem()
  {
    atomID[0] = -1;
    atomID[1] = -1;
    p[0] = NULL;
    p[1] = NULL;
  }

inline BondElem::BondElem(const Bond *a)
  {
    atomID[0] = a->atom1;
    atomID[1] = a->atom2;
    bondType = a->bond_type;
  }

inline BondElem::BondElem(AtomID atom0, AtomID atom1)
  {
    if (atom0 > atom1) {  // Swap end atoms so lowest is first!
      AtomID tmp = atom1; atom1 = atom0; atom0 = tmp; 
    }
    atomID[0] = atom0;
    atomID[1] = atom1;
  }

inline int BondElem::operator==(const BondElem &a) const
  {
    return (a.atomID[0] == atomID[0] && a.atomID[1] == atomID[1]);
  }

inline int BondElem::operator<(const BondElem &a) const
  {
    return (atomID[0] < a.atomID[0] ||
            (atomID[0] == a.atomID[0] &&
            (atomID[1] < a.atomID[1]) ));
  }

#endif

