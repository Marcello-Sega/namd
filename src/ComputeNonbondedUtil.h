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

// function arguments
#define ARGS	Position* p[2], Force* ff[2], \
		AtomProperties* a[2], \
		int numAtoms[2], BigReal *reduction

#define FARGS	Position* p[2], Force* ff[2], \
		Force* fullf [2], \
		AtomProperties* a[2], \
		int numAtoms[2], BigReal *reduction

class ComputeNonbondedUtil {

public:

  ComputeNonbondedUtil() {}
  virtual ~ComputeNonbondedUtil() {}
  static void select(void);

  static void (*calcPair)(ARGS);
  static void (*calcSelf)(ARGS);
  static void (*calcExcl)(const Position &,
		Force &, Force &,
		const AtomProperties &, const AtomProperties &,
		int, BigReal*);

  static void (*calcFullPair)(FARGS);
  static void (*calcFullSelf)(FARGS);
  static void (*calcFullExcl)(const Position &,
		Force &, Force &, Force &, Force &,
		const AtomProperties &, const AtomProperties &,
		int, BigReal*);

  enum { electEnergyIndex, fullElectEnergyIndex, vdwEnergyIndex,
	 virialIndex, fullElectVirialIndex, reductionDataSize };
  static void registerReductionData(ReductionMgr*);
  static void submitReductionData(BigReal*,ReductionMgr*,int);
  static void unregisterReductionData(ReductionMgr*);

  static Real cutoff;
  static BigReal cutoff2;
  static BigReal groupcutoff2;
  static BigReal dielectric_1;
  static const LJTable* ljTable;
  static const Molecule* mol;
  static BigReal scale14;
  static Real switchOn;
  static BigReal switchOn_1;
  static BigReal switchOn2;
  static BigReal c0;
  static BigReal c1;
  static BigReal c3;
  static BigReal c5;
  static BigReal c6;
  static BigReal d0;

  // No splitting
  static void calc_pair (ARGS);
  static void calc_pair_fullelect (FARGS);
  static void calc_self (ARGS);
  static void calc_self_fullelect (FARGS);
  static void calc_excl (
			const Position & p_ij,
			Force & f_i, Force & f_j,
			const AtomProperties & a_i, const AtomProperties & a_j,
			int m14, BigReal *reduction);
  static void calc_excl_fullelect (
			const Position & p_ij,
			Force & f_i, Force & f_j,
			Force & fullf_i, Force & fullf_j,
			const AtomProperties & a_i, const AtomProperties & a_j,
			int m14, BigReal *reduction);

  // C1 Splitting
  static void calc_pair_c1 (ARGS);
  static void calc_pair_fullelect_c1 (FARGS);
  static void calc_self_c1 (ARGS);
  static void calc_self_fullelect_c1 (FARGS);
  static void calc_excl_c1 (
			const Position & p_ij,
			Force & f_i, Force & f_j,
			const AtomProperties & a_i, const AtomProperties & a_j,
			int m14, BigReal *reduction);
  static void calc_excl_fullelect_c1 (
			const Position & p_ij,
			Force & f_i, Force & f_j,
			Force & fullf_i, Force & fullf_j,
			const AtomProperties & a_i, const AtomProperties & a_j,
			int m14, BigReal *reduction);

  // XPLOR Splitting
  static void calc_pair_xplor (ARGS);
  static void calc_pair_fullelect_xplor (FARGS);
  static void calc_self_xplor (ARGS);
  static void calc_self_fullelect_xplor (FARGS);
  static void calc_excl_xplor (
			const Position & p_ij,
			Force & f_i, Force & f_j,
			const AtomProperties & a_i, const AtomProperties & a_j,
			int m14, BigReal *reduction);
  static void calc_excl_fullelect_xplor (
			const Position & p_ij,
			Force & f_i, Force & f_j,
			Force & fullf_i, Force & fullf_j,
			const AtomProperties & a_i, const AtomProperties & a_j,
			int m14, BigReal *reduction);

};

#endif
/***************************************************************************
 * RCS INFORMATION:
 *
 *	$RCSfile: ComputeNonbondedUtil.h,v $
 *	$Author: nealk $	$Locker:  $		$State: Exp $
 *	$Revision: 1.1010 $	$Date: 1997/05/15 17:43:49 $
 *
 ***************************************************************************
 * REVISION HISTORY:
 *
 * $Log: ComputeNonbondedUtil.h,v $
 * Revision 1.1010  1997/05/15 17:43:49  nealk
 * Merged Pair and Self to use same headers.
 *
 * Revision 1.1009  1997/05/13 18:30:47  nealk
 * Removed ComputeNonbondedHack.h!
 * Reduced a lot of code in Util and Base.
 * ComputeNonbondedBase.h now only contains the function definitions.
 * The only heavy macro areas are in Util.C (determining which Base.h to define)
 * and Base.h (where the functions are defined).
 *
 * Revision 1.1008  1997/05/05 16:39:00  nealk
 * Corrected cutoff value used with hydrogen grouping.  (groupcutoff2)
 *
 * Revision 1.1007  1997/04/08 07:08:27  ari
 * Modification for dynamic loadbalancing - moving computes
 * Still bug in new computes or usage of proxies/homepatches.
 * Works if ldbStrategy is none as before.
 *
 * Revision 1.1006  1997/03/20 23:53:44  ari
 * Some changes for comments. Copyright date additions.
 * Hooks for base level update of Compute objects from ComputeMap
 * by ComputeMgr.  Useful for new compute migration functionality.
 *
 * Revision 1.1005  1997/03/16 22:56:31  jim
 * Added virial calculation for all bonded forces.
 *
 * Revision 1.1004  1997/03/14 23:18:14  jim
 * Implemented C1 splitting for long-range electrostatics.
 * Energies on sub-timesteps are incorrect.  (Aren't they always?)
 *
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

