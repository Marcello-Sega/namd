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
/***************************************************************************
 * RCS INFORMATION:
 *
 *	$RCSfile: ComputeBonds.inl,v $
 *	$Author: ari $	$Locker:  $		$State: Exp $
 *	$Revision: 1.777 $	$Date: 1997/01/17 19:35:38 $
 *
 ***************************************************************************
 * REVISION HISTORY:
 *
 * $Log: ComputeBonds.inl,v $
 * Revision 1.777  1997/01/17 19:35:38  ari
 * Internal CVS leveling release.  Start development code work
 * at 1.777.1.1.
 *
 * Revision 1.1  1997/01/14 15:29:47  nealk
 * Initial revision
 *
 * Revision 1.2  1996/12/03 15:15:40  nealk
 * Removed tons-o-debugging.
 *
 * Revision 1.1  1996/12/03 14:53:42  nealk
 * Initial revision
 *
 *
 ***************************************************************************/

