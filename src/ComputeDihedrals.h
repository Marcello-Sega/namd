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

#ifndef COMPUTEDIHEDRALS_H
#define COMPUTEDIHEDRALS_H

#include "ComputeHomeTuples.h"
class ReductionMgr;
class Molecule;

class DihedralElem {
public:
    // ComputeHomeTuples interface
    enum { size = 4 };
    AtomID atomID[size];
    int    localIndex[size];
    TuplePatchElem *p[size];
    void computeForce(BigReal*);
    // The following is evil, but the compiler chokes otherwise. (JCP)
    static void addTuplesForAtom(void*, AtomID, Molecule*);

    // Internal data
    Index dihedralType;

  enum { dihedralEnergyIndex, reductionDataSize };
  static void registerReductionData(ReductionMgr*);
  static void submitReductionData(BigReal*,ReductionMgr*,int);
  static void unregisterReductionData(ReductionMgr*);

  inline DihedralElem();
  inline DihedralElem(const Dihedral *a);
  inline DihedralElem(AtomID atom0, AtomID atom1, AtomID atom2, AtomID atom3);
  ~DihedralElem() {};

  inline int operator==(const DihedralElem &a) const;
  inline int operator<(const DihedralElem &a) const;
};

class ComputeDihedrals : public ComputeHomeTuples<DihedralElem>
{
public:

  ComputeDihedrals(ComputeID c) : ComputeHomeTuples<DihedralElem>(c) { ; }

};

#include "ComputeDihedrals.inl"

#endif
/***************************************************************************
 * RCS INFORMATION:
 *
 *	$RCSfile: ComputeDihedrals.h,v $
 *	$Author: ari $	$Locker:  $		$State: Exp $
 *	$Revision: 1.777 $	$Date: 1997/01/17 19:35:40 $
 *
 ***************************************************************************
 * REVISION HISTORY:
 *
 * $Log: ComputeDihedrals.h,v $
 * Revision 1.777  1997/01/17 19:35:40  ari
 * Internal CVS leveling release.  Start development code work
 * at 1.777.1.1.
 *
 * Revision 1.3  1997/01/16 00:55:54  jim
 * Added reduction of energies from ComputeHomeTuples objects, except
 * for ComputeNonbondedExcl which only reports 0 energy.
 * Some problems with ReductionMgr are apparent, but it still runs.
 *
 * Revision 1.2  1997/01/14 15:29:47  nealk
 * Moved "inline" functions into .inl file.
 *
 * Revision 1.1  1996/11/26 16:33:35  nealk
 * Initial revision
 *
 * Revision 1.7  1996/11/22 01:45:28  jim
 * restructured, fixed bugs, now seems to work
 *
 * Revision 1.6  1996/11/19 06:58:37  jim
 * first compiling templated version, needed ugly void* hack
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

