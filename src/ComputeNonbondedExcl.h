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

#ifndef COMPUTENONBONDEDEXCL_H
#define COMPUTENONBONDEDEXCL_H

#include "ComputeHomeTuples.h"
#include "ComputeNonbondedUtil.h"
class Molecule;

class NonbondedExclElem : public ComputeNonbondedUtil {
public:
    // ComputeHomeTuples interface
    enum { size = 2 };
    AtomID atomID[size];
    int    localIndex[size];
    TuplePatchElem *p[size];
    void computeForce(BigReal*);
    // The following is evil, but the compiler chokes otherwise. (JCP)
    static void addTuplesForAtom(void*, AtomID, Molecule*);

    // Internal data
    Index modified;

  NonbondedExclElem() {
    atomID[0] = -1;
    atomID[1] = -1;
    p[0] = NULL;
    p[1] = NULL;
  }
  NonbondedExclElem(const Exclusion *a) {
    atomID[0] = a->atom1;
    atomID[1] = a->atom2;
    modified = a->modified;
  }

  NonbondedExclElem(AtomID atom0, AtomID atom1) {
    if (atom0 > atom1) {  // Swap end atoms so lowest is first!
      AtomID tmp = atom1; atom1 = atom0; atom0 = tmp; 
    }
    atomID[0] = atom0;
    atomID[1] = atom1;
  }
  ~NonbondedExclElem() {};

  int operator==(const NonbondedExclElem &a) const {
    return (a.atomID[0] == atomID[0] && a.atomID[1] == atomID[1]);
  }

  int operator<(const NonbondedExclElem &a) const {
    return (atomID[0] < a.atomID[0] ||
            (atomID[0] == a.atomID[0] &&
            (atomID[1] < a.atomID[1]) ));
  }
};

class ComputeNonbondedExcls : public ComputeHomeTuples<NonbondedExclElem>
{
public:

  ComputeNonbondedExcls(ComputeID c) : ComputeHomeTuples<NonbondedExclElem>(c) { ; }

};

#endif
/***************************************************************************
 * RCS INFORMATION:
 *
 *	$RCSfile: ComputeNonbondedExcl.h,v $
 *	$Author: ari $	$Locker:  $		$State: Exp $
 *	$Revision: 1.777 $	$Date: 1997/01/17 19:35:54 $
 *
 ***************************************************************************
 * REVISION HISTORY:
 *
 * $Log: ComputeNonbondedExcl.h,v $
 * Revision 1.777  1997/01/17 19:35:54  ari
 * Internal CVS leveling release.  Start development code work
 * at 1.777.1.1.
 *
 * Revision 1.4  1997/01/16 20:00:03  jim
 * Added reduction calls to ComputeNonbondedSelf and ...Pair.
 * Also moved some code from ...Excl to ...Util.
 *
 * Revision 1.3  1997/01/16 00:56:01  jim
 * Added reduction of energies from ComputeHomeTuples objects, except
 * for ComputeNonbondedExcl which only reports 0 energy.
 * Some problems with ReductionMgr are apparent, but it still runs.
 *
 * Revision 1.2  1996/12/06 06:56:11  jim
 * cleaned up and renamed a bit, now it works
 *
 * Revision 1.1  1996/12/03 17:17:03  nealk
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

