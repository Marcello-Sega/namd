/**
***  Copyright (c) 1995, 1996, 1997, 1998, 1999, 2000 by
***  The Board of Trustees of the University of Illinois.
***  All rights reserved.
**/

/*
   Common operations for ComputeNonbonded classes
*/

#include "ComputeNonbondedUtil.h"
#include "SimParameters.h"
#include "Node.h"
#include "Molecule.h"
#include "LJTable.h"
#include "ReductionMgr.h"

Bool		ComputeNonbondedUtil::commOnly;
Bool		ComputeNonbondedUtil::fixedAtomsOn;
Real            ComputeNonbondedUtil::cutoff;
BigReal         ComputeNonbondedUtil::cutoff2;
BigReal         ComputeNonbondedUtil::groupcutoff2;
BigReal         ComputeNonbondedUtil::dielectric_1;
const LJTable*  ComputeNonbondedUtil::ljTable;
const Molecule* ComputeNonbondedUtil::mol;
BigReal         ComputeNonbondedUtil::scaling;
BigReal         ComputeNonbondedUtil::scale14;
Real            ComputeNonbondedUtil::switchOn;
BigReal         ComputeNonbondedUtil::switchOn_1;
BigReal         ComputeNonbondedUtil::switchOn2;
BigReal         ComputeNonbondedUtil::c0;
BigReal         ComputeNonbondedUtil::c1;
BigReal         ComputeNonbondedUtil::c3;
BigReal         ComputeNonbondedUtil::c5;
BigReal         ComputeNonbondedUtil::c6;
BigReal         ComputeNonbondedUtil::c7;
BigReal         ComputeNonbondedUtil::c8;
// BigReal         ComputeNonbondedUtil::d0;

BigReal		ComputeNonbondedUtil::ewaldcof;
BigReal		ComputeNonbondedUtil::pi_ewaldcof;

void (*ComputeNonbondedUtil::calcPair)(nonbonded *);
void (*ComputeNonbondedUtil::calcSelf)(nonbonded *);

void (*ComputeNonbondedUtil::calcFullPair)(nonbonded *);
void (*ComputeNonbondedUtil::calcFullSelf)(nonbonded *);

void (*ComputeNonbondedUtil::calcSlowPair)(nonbonded *);
void (*ComputeNonbondedUtil::calcSlowSelf)(nonbonded *);

void ComputeNonbondedUtil::submitReductionData(BigReal *data, SubmitReduction *reduction)
{
  reduction->item(REDUCTION_ELECT_ENERGY) += data[electEnergyIndex];
  reduction->item(REDUCTION_ELECT_ENERGY_SLOW) += data[fullElectEnergyIndex];
  reduction->item(REDUCTION_LJ_ENERGY) += data[vdwEnergyIndex];
  ADD_TENSOR(reduction,REDUCTION_VIRIAL_NBOND,data,virialIndex);
  ADD_TENSOR(reduction,REDUCTION_VIRIAL_SLOW,data,fullElectVirialIndex);
  reduction->item(REDUCTION_COMPUTE_CHECKSUM) += 1.;
}

#ifdef WIN32
extern "C" {
  double erfc(double) { return 0.; }
}
#else
//  Not in KCC's math.h
extern "C" {
  extern double erfc(double);
}
#endif

void ComputeNonbondedUtil::select(void)
{
  SimParameters * simParams = Node::Object()->simParameters;

  commOnly = simParams->commOnly;
  fixedAtomsOn = simParams->fixedAtomsOn;

  cutoff = simParams->cutoff;
  cutoff2 = cutoff*cutoff;

  // we add slightly more than 2 angstroms to get the same numbers.
  // don't know why...ask jim... :-)
  const BigReal &hcutoff = simParams->hgroupCutoff;
  groupcutoff2 = (cutoff+hcutoff)*(cutoff+hcutoff);

  dielectric_1 = 1.0/simParams->dielectric;
  ljTable = LJTable::Instance();
  mol = Node::Object()->molecule;
  scaling = simParams->nonbondedScaling;
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
    // d0 = 1.0/(cutoff-switchOn);
    switchOn2 = switchOn*switchOn;
    c0 = 1.0/(cutoff2-switchOn2);
  }
  else
  {
    switchOn = cutoff;
    switchOn_1 = 1.0/switchOn;
    // d0 = 0.;  // avoid division by zero
    switchOn2 = switchOn*switchOn;
    c0 = 0.;  // avoid division by zero
  }
  c1 = c0*c0*c0;
  c3 = c1 * 4.0;
  c5 = 1/cutoff2;
  c6 = -4 * c5;
  c7 = 0.5 / ( cutoff * cutoff2 );
  c8 = 1.5 / cutoff;

  int PMEOn = simParams->PMEOn;

  if ( PMEOn ) {
    ewaldcof = simParams->PMEEwaldCoefficient;
    BigReal TwoBySqrtPi = 1.12837916709551;
    pi_ewaldcof = TwoBySqrtPi * ewaldcof;
  }

  if ( ! ( simParams->fullDirectOn || simParams->FMAOn || PMEOn ) )
  {
  	calcFullPair = 0;
  	calcSlowPair = 0;
  	calcPair = calc_pair;

  	calcFullSelf = 0;
  	calcSlowSelf = 0;
  	calcSelf = calc_self;
  }
  else switch ( simParams->longSplitting )
  {
    case XPLOR:
	if ( PMEOn ) calcFullPair = calc_pair_fullelect_pme_xplor;
  	else calcFullPair = calc_pair_fullelect_xplor;
	if ( PMEOn ) calcSlowPair = calc_pair_slow_fullelect_pme_xplor;
  	else calcSlowPair = calc_pair_slow_fullelect_xplor;
  	calcPair = calc_pair_xplor;

	if ( PMEOn ) calcFullSelf = calc_self_fullelect_pme_xplor;
  	else calcFullSelf = calc_self_fullelect_xplor;
	if ( PMEOn ) calcSlowSelf = calc_self_slow_fullelect_pme_xplor;
  	else calcSlowSelf = calc_self_slow_fullelect_xplor;
  	calcSelf = calc_self_xplor;
    	break;

    case C1:
	if ( PMEOn ) calcFullPair = calc_pair_fullelect_pme_c1;
  	else calcFullPair = calc_pair_fullelect_c1;
	if ( PMEOn ) calcSlowPair = calc_pair_slow_fullelect_pme_c1;
  	else calcSlowPair = calc_pair_slow_fullelect_c1;
  	calcPair = calc_pair_c1;

	if ( PMEOn ) calcFullSelf = calc_self_fullelect_pme_c1;
  	else calcFullSelf = calc_self_fullelect_c1;
	if ( PMEOn ) calcSlowSelf = calc_self_slow_fullelect_pme_c1;
  	else calcSlowSelf = calc_self_slow_fullelect_c1;
  	calcSelf = calc_self_c1;
    	break;

    case SHARP:
    NAMD_die("Sorry, SHARP splitting not supported.");
    break;

    default:
    NAMD_die("Unknown splitting type found!");

  }
}

// clear all
// define interaction type (pair or self)
#define NBPAIR	1
#define NBSELF	2
// define electrostatics
#undef FULLELECT
#define FULLELECT_NOCORRECTION	1
#define FULLELECT_PME		2
// define splitting function
#define SPLIT_NONE	1
#define SPLIT_C1	2
#define SPLIT_XPLOR	3

// (3) BEGIN SPLITTING
#undef SPLIT_TYPE
#define SPLIT_TYPE SPLIT_NONE
//   (2) BEGIN PAIR / SELF
#undef  NBTYPE
#define NBTYPE NBPAIR
//     (1) BEGIN FULLELECT
#define SLOWONLY
#define FULLELECT FULLELECT_NOCORRECTION
#include "ComputeNonbondedBase.h"
#undef FULLELECT
#define FULLELECT FULLELECT_PME
#include "ComputeNonbondedBase.h"
#undef FULLELECT
#undef SLOWONLY
#define FULLELECT FULLELECT_NOCORRECTION
#include "ComputeNonbondedBase.h"
#undef FULLELECT
#define FULLELECT FULLELECT_PME
#include "ComputeNonbondedBase.h"
#undef FULLELECT
#include "ComputeNonbondedBase.h"
//     (1) END FULLELECT
#undef  NBTYPE
#define NBTYPE NBSELF
//     (1) BEGIN FULLELECT
#define SLOWONLY
#define FULLELECT FULLELECT_NOCORRECTION
#include "ComputeNonbondedBase.h"
#undef FULLELECT
#define FULLELECT FULLELECT_PME
#include "ComputeNonbondedBase.h"
#undef FULLELECT
#undef SLOWONLY
#define FULLELECT FULLELECT_NOCORRECTION
#include "ComputeNonbondedBase.h"
#undef FULLELECT
#define FULLELECT FULLELECT_PME
#include "ComputeNonbondedBase.h"
#undef FULLELECT
#include "ComputeNonbondedBase.h"
//     (1) END FULLELECT
//   (2) END PAIR / SELF

#undef SPLIT_TYPE
#define SPLIT_TYPE SPLIT_XPLOR
//   (2) BEGIN PAIR / SELF
#undef  NBTYPE
#define NBTYPE NBPAIR
//     (1) BEGIN FULLELECT
#define SLOWONLY
#define FULLELECT FULLELECT_NOCORRECTION
#include "ComputeNonbondedBase.h"
#undef FULLELECT
#define FULLELECT FULLELECT_PME
#include "ComputeNonbondedBase.h"
#undef FULLELECT
#undef SLOWONLY
#define FULLELECT FULLELECT_NOCORRECTION
#include "ComputeNonbondedBase.h"
#undef FULLELECT
#define FULLELECT FULLELECT_PME
#include "ComputeNonbondedBase.h"
#undef FULLELECT
#include "ComputeNonbondedBase.h"
//     (1) END FULLELECT
#undef  NBTYPE
#define NBTYPE NBSELF
//     (1) BEGIN FULLELECT
#define SLOWONLY
#define FULLELECT FULLELECT_NOCORRECTION
#include "ComputeNonbondedBase.h"
#undef FULLELECT
#define FULLELECT FULLELECT_PME
#include "ComputeNonbondedBase.h"
#undef FULLELECT
#undef SLOWONLY
#define FULLELECT FULLELECT_NOCORRECTION
#include "ComputeNonbondedBase.h"
#undef FULLELECT
#define FULLELECT FULLELECT_PME
#include "ComputeNonbondedBase.h"
#undef FULLELECT
#include "ComputeNonbondedBase.h"
//     (1) END FULLELECT
//   (2) END PAIR / SELF

#undef SPLIT_TYPE
#define SPLIT_TYPE SPLIT_C1
//   (2) BEGIN PAIR / SELF
#undef  NBTYPE
#define NBTYPE NBPAIR
//     (1) BEGIN FULLELECT
#define SLOWONLY
#define FULLELECT FULLELECT_NOCORRECTION
#include "ComputeNonbondedBase.h"
#undef FULLELECT
#define FULLELECT FULLELECT_PME
#include "ComputeNonbondedBase.h"
#undef FULLELECT
#undef SLOWONLY
#define FULLELECT FULLELECT_NOCORRECTION
#include "ComputeNonbondedBase.h"
#undef FULLELECT
#define FULLELECT FULLELECT_PME
#include "ComputeNonbondedBase.h"
#undef FULLELECT
#include "ComputeNonbondedBase.h"
//     (1) END FULLELECT
#undef  NBTYPE
#define NBTYPE NBSELF
//     (1) BEGIN FULLELECT
#define SLOWONLY
#define FULLELECT FULLELECT_NOCORRECTION
#include "ComputeNonbondedBase.h"
#undef FULLELECT
#define FULLELECT FULLELECT_PME
#include "ComputeNonbondedBase.h"
#undef FULLELECT
#undef SLOWONLY
#define FULLELECT FULLELECT_NOCORRECTION
#include "ComputeNonbondedBase.h"
#undef FULLELECT
#define FULLELECT FULLELECT_PME
#include "ComputeNonbondedBase.h"
#undef FULLELECT
#include "ComputeNonbondedBase.h"
//     (1) END FULLELECT
//   (2) END PAIR / SELF
// (3) END SPLITTING

