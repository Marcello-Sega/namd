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

#ifndef COMPUTEBOND_H
#define COMPUTEBOND_H

#include "ComputeHomeTuples.h"
class ReductionMgr;
class Molecule;

class BondElem {
public:
    // ComputeHomeTuples interface
    enum { size = 2 };
    AtomID atomID[size];
    int    localIndex[size];
    TuplePatchElem *p[size];
    void computeForce(BigReal*);
    // The following is evil, but the compiler chokes otherwise. (JCP)
    static void loadTuplesForAtom(void*, AtomID, Molecule*);

    // Internal data
    Index bondType;

  enum { bondEnergyIndex, reductionDataSize };
  static void registerReductionData(ReductionMgr*);
  static void submitReductionData(BigReal*,ReductionMgr*,int);
  static void unregisterReductionData(ReductionMgr*);

  inline BondElem();
  inline BondElem(const Bond *a);
  inline BondElem(AtomID atom0, AtomID atom1);
  ~BondElem() {};

  inline int operator==(const BondElem &a) const;
  inline int operator<(const BondElem &a) const;
};

class ComputeBonds : public ComputeHomeTuples<BondElem>
{
public:

  ComputeBonds(ComputeID c) : ComputeHomeTuples<BondElem>(c) { ; }

};

#include "ComputeBonds.inl"

#endif
/***************************************************************************
 * RCS INFORMATION:
 *
 *	$RCSfile: ComputeBonds.h,v $
 *	$Author: ari $	$Locker:  $		$State: Exp $
 *	$Revision: 1.1000 $	$Date: 1997/02/06 15:57:47 $
 *
 ***************************************************************************
 * REVISION HISTORY:
 *
 * $Log: ComputeBonds.h,v $
 * Revision 1.1000  1997/02/06 15:57:47  ari
 * Resetting CVS to merge branches back into the main trunk.
 * We will stick to main trunk development as suggested by CVS manual.
 * We will set up tags to track fixed points of development/release
 * as suggested by CVS manual - all praise the CVS manual.
 *
 * Revision 1.779  1997/02/06 15:52:51  ari
 * Updating Revision Line, getting rid of branches
 *
 * Revision 1.778.2.1  1997/02/05 22:17:58  ari
 * Added migration code - Currently the framework is
 * there with compiling code.  This version does
 * crash shortly after migration is complete.
 * Migration appears to complete, but Patches do
 * not appear to be left in a correct state.
 *
 * Revision 1.778  1997/01/28 00:30:04  ari
 * internal release uplevel to 1.778
 *
 * Revision 1.777.2.1  1997/01/27 22:44:57  ari
 * Basic Atom Migration Code added.
 * Added correct magic first line to .h files for xemacs to go to C++ mode.
 * Compiles and runs without migration turned on.
 *
 * Revision 1.777  1997/01/17 19:35:37  ari
 * Internal CVS leveling release.  Start development code work
 * at 1.777.1.1.
 *
 * Revision 1.4  1997/01/16 00:55:51  jim
 * Added reduction of energies from ComputeHomeTuples objects, except
 * for ComputeNonbondedExcl which only reports 0 energy.
 * Some problems with ReductionMgr are apparent, but it still runs.
 *
 * Revision 1.3  1997/01/14 15:29:47  nealk
 * Moved "inline" functions into .inl file.
 *
 * Revision 1.2  1996/12/03 15:15:40  nealk
 * Removed tons-o-debugging.
 *
 * Revision 1.1  1996/12/03 14:53:42  nealk
 * Initial revision
 *
 *
 ***************************************************************************/

