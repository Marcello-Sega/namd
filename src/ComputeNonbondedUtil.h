//-*-c++-*-
/***************************************************************************/
/*                                                                         */
/*              (C) Copyright 1996 The Board of Trustees of the            */
/*                          University of Illinois                         */
/*                           All Rights Reserved                           */
/*									   */
/***************************************************************************/

/***************************************************************************
 * DESCRIPTION: Common operations for ComputeNonbonded classes
 *
 ***************************************************************************/

#ifndef COMPUTENONBONDEDUTIL_H
#define COMPUTENONBONDEDUTIL_H

#include "NamdTypes.h"
class LJTable;
class ReductionMgr;
class Molecule;

class ComputeNonbondedUtil {

public:

  static void select(void);

  static void (*calcPair)(Position*[2],Force*[2],AtomProperties*[2],int[2],BigReal*);
  static void (*calcSelf)(Position*,Force*,AtomProperties*,int,BigReal*);
  static void (*calcExcl)(const Position &,
		Force &, Force &,
		const AtomProperties &, const AtomProperties &,
		int, BigReal*);

  static void (*calcFullPair)(Position*[2],Force*[2],Force*[2],AtomProperties*[2],int[2],BigReal*);
  static void (*calcFullSelf)(Position*,Force*,Force*,AtomProperties*,int,BigReal*);
  static void (*calcFullExcl)(const Position &,
		Force &, Force &, Force &, Force &,
		const AtomProperties &, const AtomProperties &,
		int, BigReal*);

  enum { electEnergyIndex, fullElectEnergyIndex, vdwEnergyIndex, reductionDataSize };
  static void registerReductionData(ReductionMgr*);
  static void submitReductionData(BigReal*,ReductionMgr*,int);
  static void unregisterReductionData(ReductionMgr*);

  static Real cutoff;
  static BigReal cutoff2;
  static BigReal dielectric_1;
  static const LJTable* ljTable;
  static const Molecule* mol;
  static BigReal scale14;
  static Real switchOn;
  static BigReal switchOn2;
  static BigReal c0;
  static BigReal c1;
  static BigReal c3;
  static BigReal c5;
  static BigReal c6;

#define DECLARATION
#undef DEFINITON

// (3) BEGIN SPLITTING
#define NOSPLIT
#undef SPLIT_XPLOR
#undef SPLIT_C1
//   (2) BEGIN PAIR / SELF / EXCL
#define NBPAIR
#undef NBSELF
#undef NBEXCL
//     (1) BEGIN FULLELECT
#define FULLELECT
#include "ComputeNonbondedBase.h"
#undef FULLELECT
#include "ComputeNonbondedBase.h"
//     (1) END FULLELECT
#undef NBPAIR
#define NBSELF
#undef NBEXCL
//     (1) BEGIN FULLELECT
#define FULLELECT
#include "ComputeNonbondedBase.h"
#undef FULLELECT
#include "ComputeNonbondedBase.h"
//     (1) END FULLELECT
#undef NBPAIR
#undef NBSELF
#define NBEXCL
//     (1) BEGIN FULLELECT
#define FULLELECT
#include "ComputeNonbondedBase.h"
#undef FULLELECT
#include "ComputeNonbondedBase.h"
//     (1) END FULLELECT
//   (2) END PAIR / SELF / EXCL

#undef NOSPLIT
#define SPLIT_XPLOR
#undef SPLIT_C1
//   (2) BEGIN PAIR / SELF / EXCL
#define NBPAIR
#undef NBSELF
#undef NBEXCL
//     (1) BEGIN FULLELECT
#define FULLELECT
#include "ComputeNonbondedBase.h"
#undef FULLELECT
#include "ComputeNonbondedBase.h"
//     (1) END FULLELECT
#undef NBPAIR
#define NBSELF
#undef NBEXCL
//     (1) BEGIN FULLELECT
#define FULLELECT
#include "ComputeNonbondedBase.h"
#undef FULLELECT
#include "ComputeNonbondedBase.h"
//     (1) END FULLELECT
#undef NBPAIR
#undef NBSELF
#define NBEXCL
//     (1) BEGIN FULLELECT
#define FULLELECT
#include "ComputeNonbondedBase.h"
#undef FULLELECT
#include "ComputeNonbondedBase.h"
//     (1) END FULLELECT
//   (2) END PAIR / SELF / EXCL

#undef NOSPLIT
#undef SPLIT_XPLOR
#define SPLIT_C1
//   (2) BEGIN PAIR / SELF / EXCL
#define NBPAIR
#undef NBSELF
#undef NBEXCL
//     (1) BEGIN FULLELECT
#define FULLELECT
#include "ComputeNonbondedBase.h"
#undef FULLELECT
#include "ComputeNonbondedBase.h"
//     (1) END FULLELECT
#undef NBPAIR
#define NBSELF
#undef NBEXCL
//     (1) BEGIN FULLELECT
#define FULLELECT
#include "ComputeNonbondedBase.h"
#undef FULLELECT
#include "ComputeNonbondedBase.h"
//     (1) END FULLELECT
#undef NBPAIR
#undef NBSELF
#define NBEXCL
//     (1) BEGIN FULLELECT
#define FULLELECT
#include "ComputeNonbondedBase.h"
#undef FULLELECT
#include "ComputeNonbondedBase.h"
//     (1) END FULLELECT
//   (2) END PAIR / SELF / EXCL
// (3) END SPLITTING

};

#endif
/***************************************************************************
 * RCS INFORMATION:
 *
 *	$RCSfile: ComputeNonbondedUtil.h,v $
 *	$Author: jim $	$Locker:  $		$State: Exp $
 *	$Revision: 1.1003 $	$Date: 1997/03/14 06:44:58 $
 *
 ***************************************************************************
 * REVISION HISTORY:
 *
 * $Log: ComputeNonbondedUtil.h,v $
 * Revision 1.1003  1997/03/14 06:44:58  jim
 * First working versions of full electrostatics splitting functions.
 *
 * Revision 1.1002  1997/02/28 04:47:06  jim
 * Full electrostatics now works with fulldirect on one node.
 *
 * Revision 1.1001  1997/02/21 20:45:13  jim
 * Eliminated multiple function for switching and modified 1-4 interactions.
 * Now assumes a switching function, but parameters are such that nothing
 * happens, same for modified 1-4.  Slight penalty for rare simulations
 * in which these features are not used, but otherwise no loss and
 * simplifies code.
 *
 * Revision 1.1000  1997/02/06 15:58:14  ari
 * Resetting CVS to merge branches back into the main trunk.
 * We will stick to main trunk development as suggested by CVS manual.
 * We will set up tags to track fixed points of development/release
 * as suggested by CVS manual - all praise the CVS manual.
 *
 * Revision 1.778  1997/01/28 00:30:26  ari
 * internal release uplevel to 1.778
 *
 * Revision 1.777.2.1  1997/01/27 22:45:08  ari
 * Basic Atom Migration Code added.
 * Added correct magic first line to .h files for xemacs to go to C++ mode.
 * Compiles and runs without migration turned on.
 *
 * Revision 1.777  1997/01/17 19:36:00  ari
 * Internal CVS leveling release.  Start development code work
 * at 1.777.1.1.
 *
 * Revision 1.7  1997/01/16 20:00:25  jim
 * Added reduction calls to ComputeNonbondedSelf and ...Pair.
 * Also moved some code from ...Excl to ...Util.
 *
 * Revision 1.6  1996/12/04 17:16:32  jim
 * ComputeNonbondedUtil::select() now caches simulation parameters
 *
 * Revision 1.5  1996/12/03 21:05:09  jim
 * added support for exclusion correction computes
 *
 * Revision 1.4  1996/11/21 01:00:15  jim
 * made methods public, got rid of friends
 *
 * Revision 1.3  1996/11/21 00:00:40  jim
 * added select(), calcPair, and calcSelf
 *
 * Revision 1.2  1996/11/20 23:16:39  jim
 * first compiling version of generic nonbonded function
 *
 * Revision 1.1  1996/10/31 22:35:04  jim
 * Initial revision
 *
 *
 ***************************************************************************/

