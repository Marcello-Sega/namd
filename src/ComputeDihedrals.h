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

#ifndef COMPUTEDIHEDRALS_H
#define COMPUTEDIHEDRALS_H

#include "ComputeHomeTuples.h"
#include "ReductionMgr.h"

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
    static void loadTuplesForAtom(void*, AtomID, Molecule*);
    static void getMoleculePointers(Molecule*, int*, int***, Dihedral**);

    // Internal data
    Index dihedralType;

  int hash() const { 
    return 0x7FFFFFFF &((atomID[0]<<24) + (atomID[1]<<16) + (atomID[2]<<8) + atomID[3]);
  }

  enum { dihedralEnergyIndex, virialXIndex, virialYIndex, virialZIndex, reductionDataSize };
  enum { reductionChecksumLabel = REDUCTION_DIHEDRAL_CHECKSUM };
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

class ComputeDihedrals : public ComputeHomeTuples<DihedralElem,Dihedral>
{
public:

  ComputeDihedrals(ComputeID c) : ComputeHomeTuples<DihedralElem,Dihedral>(c) { ; }

};

#include "ComputeDihedrals.inl"

#endif
/***************************************************************************
 * RCS INFORMATION:
 *
 *	$RCSfile: ComputeDihedrals.h,v $
 *	$Author: jim $	$Locker:  $		$State: Exp $
 *	$Revision: 1.1005 $	$Date: 1999/01/06 00:56:21 $
 *
 ***************************************************************************
 * REVISION HISTORY:
 *
 * $Log: ComputeDihedrals.h,v $
 * Revision 1.1005  1999/01/06 00:56:21  jim
 * All compute objects except DPMTA now return diagonal of virial tensor.
 *
 * Revision 1.1004  1998/11/01 23:25:45  jim
 * Added basic correctness checking: atom counts, etc.
 *
 * Revision 1.1003  1997/10/17 17:16:45  jim
 * Switched from hash tables to checklists, eliminated special exclusion code.
 *
 * Revision 1.1002  1997/03/16 22:56:25  jim
 * Added virial calculation for all bonded forces.
 *
 * Revision 1.1001  1997/03/10 17:40:04  ari
 * UniqueSet changes - some more commenting and cleanup
 *
 * Revision 1.1000  1997/02/06 15:57:51  ari
 * Resetting CVS to merge branches back into the main trunk.
 * We will stick to main trunk development as suggested by CVS manual.
 * We will set up tags to track fixed points of development/release
 * as suggested by CVS manual - all praise the CVS manual.
 *
 * Revision 1.779  1997/02/06 15:52:53  ari
 * Updating Revision Line, getting rid of branches
 *
 * Revision 1.778.2.1  1997/02/05 22:18:00  ari
 * Added migration code - Currently the framework is
 * there with compiling code.  This version does
 * crash shortly after migration is complete.
 * Migration appears to complete, but Patches do
 * not appear to be left in a correct state.
 *
 * Revision 1.778  1997/01/28 00:30:07  ari
 * internal release uplevel to 1.778
 *
 * Revision 1.777.2.1  1997/01/27 22:44:58  ari
 * Basic Atom Migration Code added.
 * Added correct magic first line to .h files for xemacs to go to C++ mode.
 * Compiles and runs without migration turned on.
 *
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

