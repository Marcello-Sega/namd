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

#include "ComputeNonbondedUtil.h"
#include "SimParameters.h"
#include "Node.h"
#include "Molecule.h"
#include "LJTable.h"
#include "ReductionMgr.h"

Real            ComputeNonbondedUtil::cutoff;
BigReal         ComputeNonbondedUtil::cutoff2;
BigReal         ComputeNonbondedUtil::groupcutoff2;
BigReal         ComputeNonbondedUtil::dielectric_1;
const LJTable*  ComputeNonbondedUtil::ljTable;
const Molecule* ComputeNonbondedUtil::mol;
BigReal         ComputeNonbondedUtil::scale14;
Real            ComputeNonbondedUtil::switchOn;
BigReal         ComputeNonbondedUtil::switchOn_1;
BigReal         ComputeNonbondedUtil::switchOn2;
BigReal         ComputeNonbondedUtil::c0;
BigReal         ComputeNonbondedUtil::c1;
BigReal         ComputeNonbondedUtil::c3;
BigReal         ComputeNonbondedUtil::c5;
BigReal         ComputeNonbondedUtil::c6;
BigReal         ComputeNonbondedUtil::d0;


void ComputeNonbondedUtil::registerReductionData(ReductionMgr *reduction)
{
  reduction->Register(REDUCTION_ELECT_ENERGY);
  reduction->Register(REDUCTION_LJ_ENERGY);
  reduction->Register(REDUCTION_VIRIAL);
}

void ComputeNonbondedUtil::submitReductionData(BigReal *data, ReductionMgr *reduction, int seq)
{
  reduction->submit(seq, REDUCTION_ELECT_ENERGY, data[electEnergyIndex]
					+ data[fullElectEnergyIndex]);
  reduction->submit(seq, REDUCTION_LJ_ENERGY, data[vdwEnergyIndex]);
  reduction->submit(seq, REDUCTION_VIRIAL, data[virialIndex]
					+ data[fullElectVirialIndex]);
}

void ComputeNonbondedUtil::unregisterReductionData(ReductionMgr *reduction)
{
  reduction->unRegister(REDUCTION_ELECT_ENERGY);
  reduction->unRegister(REDUCTION_LJ_ENERGY);
  reduction->unRegister(REDUCTION_VIRIAL);
}


void ComputeNonbondedUtil::select(void)
{
  SimParameters * simParams = Node::Object()->simParameters;

  cutoff = simParams->cutoff;
  cutoff2 = cutoff*cutoff;
  // we add slightly more than 2 angstroms to get the same numbers.
  // don't know why...ask jim... :-)
  groupcutoff2 = (cutoff+2.5)*(cutoff+2.5);
  dielectric_1 = 1.0/simParams->dielectric;
  ljTable = LJTable::Instance();
  mol = Node::Object()->molecule;
  if ( simParams->exclude == SCALED14 )
  {
    scale14 = simParams->scale14;
  }
  else
  {
    scale14 = 1.;
  }
  if ( simParams->switchingActive )
  {
    switchOn = simParams->switchingDist;
    switchOn_1 = 1.0/switchOn;
    d0 = 1.0/(cutoff-switchOn);
    switchOn2 = switchOn*switchOn;
    c0 = 1.0/(cutoff2-switchOn2);
  }
  else
  {
    switchOn = cutoff;
    switchOn_1 = 1.0/switchOn;
    d0 = 0.;  // avoid division by zero
    switchOn2 = switchOn*switchOn;
    c0 = 0.;  // avoid division by zero
  }
  c1 = c0*c0*c0;
  c3 = c1 * 4.0;
  c5 = 1/cutoff2;
  c6 = -4 * c5;

  if ( ! ( simParams->fullDirectOn || simParams->FMAOn ) )
  {
  	calcFullPair = 0;
  	calcPair = calc_pair;

  	calcFullSelf = 0;
  	calcSelf = calc_self;

  	calcFullExcl = 0;
  	calcExcl = calc_excl;
  }
  else switch ( simParams->longSplitting )
  {
    case XPLOR:
  	calcFullPair = calc_pair_fullelect_xplor;
  	calcPair = calc_pair_xplor;

  	calcFullSelf = calc_self_fullelect_xplor;
  	calcSelf = calc_self_xplor;

  	calcFullExcl = calc_excl_fullelect_xplor;
  	calcExcl = calc_excl_xplor;
    	break;

    case C1:
  	calcFullPair = calc_pair_fullelect_c1;
  	calcPair = calc_pair_c1;

  	calcFullSelf = calc_self_fullelect_c1;
  	calcSelf = calc_self_c1;

  	calcFullExcl = calc_excl_fullelect_c1;
  	calcExcl = calc_excl_c1;
	break;

    case SKEEL:
    NAMD_die("Sorry, SKEEL splitting not supported.");
    break;

    case SHARP:
    NAMD_die("Sorry, SHARP splitting not supported.");
    break;

    default:
    NAMD_die("Unknown splitting type found!");

  }
}

#define NBPAIR
#undef NBSELF
#undef NBEXCL

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


/***************************************************************************
 * RCS INFORMATION:
 *
 *	$RCSfile: ComputeNonbondedUtil.C,v $
 *	$Author: nealk $	$Locker:  $		$State: Exp $
 *	$Revision: 1.1011 $	$Date: 1997/05/15 17:43:48 $
 *
 ***************************************************************************
 * REVISION HISTORY:
 *
 * $Log: ComputeNonbondedUtil.C,v $
 * Revision 1.1011  1997/05/15 17:43:48  nealk
 * Merged Pair and Self to use same headers.
 *
 * Revision 1.1010  1997/05/13 18:30:47  nealk
 * Removed ComputeNonbondedHack.h!
 * Reduced a lot of code in Util and Base.
 * ComputeNonbondedBase.h now only contains the function definitions.
 * The only heavy macro areas are in Util.C (determining which Base.h to define)
 * and Base.h (where the functions are defined).
 *
 * Revision 1.1009  1997/05/09 18:24:23  nealk
 * 1. Added hydrogen grouping code to improve performance in ComputeNonbondedBase
 *    CODE ONLY WORKS WITH HYDROGEN GROUPING!
 * 2. Increased the hydrogen group cutoff side from 2A to 2.5A -- 2A gave
 *    fractionally different values after 100 iterations.  2.5A gives same numbers.
 * 3. Made migration by hydrogen grouping the default in SimParameters.
 *
 * Revision 1.1008  1997/05/05 16:38:59  nealk
 * Corrected cutoff value used with hydrogen grouping.  (groupcutoff2)
 *
 * Revision 1.1007  1997/04/04 23:34:18  milind
 * Got NAMD2 to run on Origin2000.
 * Included definitions of class static variables in C files.
 * Fixed alignment bugs by using memcpy instead of assignment in
 * pack and unpack.
 *
 * Revision 1.1006  1997/03/20 23:53:43  ari
 * Some changes for comments. Copyright date additions.
 * Hooks for base level update of Compute objects from ComputeMap
 * by ComputeMgr.  Useful for new compute migration functionality.
 *
 * Revision 1.1005  1997/03/16 22:56:29  jim
 * Added virial calculation for all bonded forces.
 *
 * Revision 1.1004  1997/03/14 23:18:13  jim
 * Implemented C1 splitting for long-range electrostatics.
 * Energies on sub-timesteps are incorrect.  (Aren't they always?)
 *
 * Revision 1.1003  1997/03/14 06:44:56  jim
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
 * Revision 1.1000  1997/02/06 15:58:13  ari
 * Resetting CVS to merge branches back into the main trunk.
 * We will stick to main trunk development as suggested by CVS manual.
 * We will set up tags to track fixed points of development/release
 * as suggested by CVS manual - all praise the CVS manual.
 *
 * Revision 1.778  1997/01/28 00:30:25  ari
 * internal release uplevel to 1.778
 *
 * Revision 1.777  1997/01/17 19:35:59  ari
 * Internal CVS leveling release.  Start development code work
 * at 1.777.1.1.
 *
 * Revision 1.6  1997/01/16 20:00:23  jim
 * Added reduction calls to ComputeNonbondedSelf and ...Pair.
 * Also moved some code from ...Excl to ...Util.
 *
 * Revision 1.5  1996/12/04 17:16:32  jim
 * ComputeNonbondedUtil::select() now caches simulation parameters
 *
 * Revision 1.4  1996/12/03 21:05:09  jim
 * added support for exclusion correction computes
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

