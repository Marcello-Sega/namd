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

#include "ComputeHomeTuples.h"
class ReductionMgr;
class Molecule;

class AngleElem {
public:
    // ComputeHomeTuples interface
    enum { size = 3 };
    AtomID atomID[size];
    int    localIndex[size];
    TuplePatchElem *p[size];
    void computeForce(BigReal*);
    // The following is evil, but the compiler chokes otherwise. (JCP)
    static void addTuplesForAtom(void*, AtomID, Molecule*);

    // Internal data
    Index angleType;

  enum { angleEnergyIndex, reductionDataSize };
  static void registerReductionData(ReductionMgr*);
  static void submitReductionData(BigReal*,ReductionMgr*,int);
  static void unregisterReductionData(ReductionMgr*);

  inline AngleElem();
  inline AngleElem(const Angle *a);
  inline AngleElem(AtomID atom0, AtomID atom1, AtomID atom2);
  ~AngleElem() {};

  inline int operator==(const AngleElem &a) const;
  inline int operator<(const AngleElem &a) const;
};

class ComputeAngles : public ComputeHomeTuples<AngleElem>
{
public:

  ComputeAngles(ComputeID c) : ComputeHomeTuples<AngleElem>(c) { ; }

};

#include "ComputeAngles.inl"

#endif
/***************************************************************************
 * RCS INFORMATION:
 *
 *	$RCSfile: ComputeAngles.h,v $
 *	$Author: ari $	$Locker:  $		$State: Exp $
 *	$Revision: 1.777 $	$Date: 1997/01/17 19:35:35 $
 *
 ***************************************************************************
 * REVISION HISTORY:
 *
 * $Log: ComputeAngles.h,v $
 * Revision 1.777  1997/01/17 19:35:35  ari
 * Internal CVS leveling release.  Start development code work
 * at 1.777.1.1.
 *
 * Revision 1.9  1997/01/16 00:55:49  jim
 * Added reduction of energies from ComputeHomeTuples objects, except
 * for ComputeNonbondedExcl which only reports 0 energy.
 * Some problems with ReductionMgr are apparent, but it still runs.
 *
 * Revision 1.8  1997/01/14 15:29:47  nealk
 * Moved "include" functions into .inl file.
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

