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
#include "InfoStream.h"
#include <stdio.h>

Bool		ComputeNonbondedUtil::commOnly;
Bool		ComputeNonbondedUtil::fixedAtomsOn;
Real            ComputeNonbondedUtil::cutoff;
BigReal         ComputeNonbondedUtil::cutoff2;
BigReal         ComputeNonbondedUtil::groupcutoff2;
BigReal         ComputeNonbondedUtil::dielectric_1;
const LJTable*  ComputeNonbondedUtil::ljTable = 0;
const Molecule* ComputeNonbondedUtil::mol;
BigReal		ComputeNonbondedUtil::r2_delta;
BigReal		ComputeNonbondedUtil::r2_delta_1;
BigReal*	ComputeNonbondedUtil::table_alloc = 0;
BigReal*	ComputeNonbondedUtil::fast_table;
BigReal*	ComputeNonbondedUtil::scor_table;
BigReal*	ComputeNonbondedUtil::slow_table;
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

void (*ComputeNonbondedUtil::calcPair)(nonbonded *) = calc_pair;
void (*ComputeNonbondedUtil::calcSelf)(nonbonded *) = calc_self;

void (*ComputeNonbondedUtil::calcFullPair)(nonbonded *) = calc_pair_fullelect;
void (*ComputeNonbondedUtil::calcFullSelf)(nonbonded *) = calc_self_fullelect;

void (*ComputeNonbondedUtil::calcSlowPair)(nonbonded *) = calc_pair_slow_fullelect;
void (*ComputeNonbondedUtil::calcSlowSelf)(nonbonded *) = calc_self_slow_fullelect;

// define splitting function
#define SPLIT_NONE	1
#define SPLIT_SHIFT	2
#define SPLIT_C1	3
#define SPLIT_XPLOR	4

void ComputeNonbondedUtil::submitReductionData(BigReal *data, SubmitReduction *reduction)
{
  reduction->item(REDUCTION_EXCLUSION_CHECKSUM) += data[exclChecksumIndex];
  reduction->item(REDUCTION_ELECT_ENERGY) += data[electEnergyIndex];
  reduction->item(REDUCTION_ELECT_ENERGY_SLOW) += data[fullElectEnergyIndex];
  reduction->item(REDUCTION_LJ_ENERGY) += data[vdwEnergyIndex];
  ADD_TENSOR(reduction,REDUCTION_VIRIAL_NBOND,data,virialIndex);
  ADD_TENSOR(reduction,REDUCTION_VIRIAL_SLOW,data,fullElectVirialIndex);
  reduction->item(REDUCTION_COMPUTE_CHECKSUM) += 1.;
}

#ifndef _IA64
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
  if ( ! ljTable ) ljTable = new LJTable;
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

  const int PMEOn = simParams->PMEOn;

  if ( PMEOn ) {
    ewaldcof = simParams->PMEEwaldCoefficient;
    BigReal TwoBySqrtPi = 1.12837916709551;
    pi_ewaldcof = TwoBySqrtPi * ewaldcof;
  }

  int splitType = SPLIT_NONE;
  if ( simParams->switchingActive ) splitType = SPLIT_SHIFT;
  if ( simParams->fullDirectOn || simParams->FMAOn || PMEOn ) {
    switch ( simParams->longSplitting ) {
      case C1:
      splitType = SPLIT_C1;
      break;

      case XPLOR:
      NAMD_die("Sorry, XPLOR splitting not supported.");
      break;

      case SHARP:
      NAMD_die("Sorry, SHARP splitting not supported.");
      break;

      default:
      NAMD_die("Unknown splitting type found!");

    }
  }

  BigReal r2_tol = 0.1;
  
  r2_delta = 1.0;
  while ( r2_delta > r2_tol ) r2_delta /= 2.0;
  r2_delta_1 = 1.0 / r2_delta;

  if ( ! CkMyPe() ) {
    iout << iINFO << "COULOMB TABLE R-SQUARED SPACING: " <<
				r2_delta << "\n" << endi;
  }

  int i;
  int n = (int)(cutoff2 / r2_delta) + 4;
  n += ( (n % 2) ? 0 : 1 );  // ensure that n is odd

  if ( ! CkMyPe() ) {
    iout << iINFO << "COULOMB TABLE SIZE: " <<
				n << " POINTS\n" << endi;
  }

  if ( table_alloc ) delete [] table_alloc;
  table_alloc = new BigReal[12*n+40];
  BigReal *table_align = table_alloc;
  while ( ((long)table_align) % 32 ) ++table_align;
  fast_table = table_align;
  scor_table = table_align + 4*n;
  slow_table = table_align + 8*n;
  BigReal *fast_i = fast_table + 4;
  BigReal *scor_i = scor_table + 4;
  BigReal *slow_i = slow_table + 4;

  // fill in the rest of the table
  for ( i=1; i<n; ++i ) {

    const BigReal r2 = r2_delta * i;

    const BigReal r = sqrt(r2);
    const BigReal r_1 = 1.0/r;

    // fast_ is defined as (full_ - slow_)
    // corr_ and fast_ are both zero at the cutoff, full_ is not
    // all three are approx 1/r at short distances

    // for actual interpolation, we use fast_ for fast forces and
    // scor_ = slow_ + corr_ - full_ and slow_ for slow forces
    // since these last two are of small magnitude

    BigReal fast_energy, fast_gradient;
    BigReal scor_energy, scor_gradient;
    BigReal slow_energy, slow_gradient;

    // corr_ is PME direct sum, or similar correction term
    // corr_energy is multiplied by r until later
    // corr_gradient is multiplied by -r^2 until later
    BigReal corr_energy, corr_gradient;

    if ( PMEOn ) {
      BigReal tmp_a = r * ewaldcof;
      BigReal tmp_b = erfc(tmp_a);
      corr_energy = tmp_b;
      corr_gradient = pi_ewaldcof*exp(-(tmp_a*tmp_a))*r + tmp_b;
    } else {
      corr_energy = corr_gradient = 0;
    }

    switch(splitType) {
      case SPLIT_NONE:
        fast_energy = 1.0/r;
        fast_gradient = -1.0/r2;
        scor_energy = scor_gradient = 0;
        slow_energy = slow_gradient = 0;
	break;
      case SPLIT_SHIFT: {
	BigReal shiftVal = r2/cutoff2 - 1.0;
	shiftVal *= shiftVal;
	BigReal dShiftVal = 2.0 * (r2/cutoff2 - 1.0) * 2.0*r/cutoff2;
        fast_energy = shiftVal/r;
        fast_gradient = dShiftVal/r - shiftVal/r2;
        scor_energy = scor_gradient = 0;
        slow_energy = slow_gradient = 0;
        } 
	break;
      case SPLIT_C1:
	// calculate actual energy and gradient
	slow_energy = 0.5/cutoff * (3.0 - (r2/cutoff2));
	slow_gradient = -1.0/cutoff2 * (r/cutoff);
	// calculate scor from slow and corr
	scor_energy = slow_energy + (corr_energy - 1.0)/r;
	scor_gradient = slow_gradient - (corr_gradient - 1.0)/r2;
	// calculate fast from slow
	fast_energy = 1.0/r - slow_energy;
	fast_gradient = -1.0/r2 - slow_gradient;
	break;
    }

    // foo_gradient is calculated as ( d foo_energy / d r )
    // and now divided by 2r to get ( d foo_energy / d r2 )

    fast_gradient *= 0.5 * r_1;
    scor_gradient *= 0.5 * r_1;
    slow_gradient *= 0.5 * r_1;

    // let modf be 1 if excluded, 1-scale14 if modified, 0 otherwise,
    // add scor_ - modf * slow_ to slow terms and
    // add fast_ - modf * fast_ to fast terms.

    *(fast_i++) = fast_energy;
    *(fast_i++) = fast_gradient;
    *(fast_i++) = 0;
    *(fast_i++) = 0;
    *(scor_i++) = scor_energy;
    *(scor_i++) = scor_gradient;
    *(scor_i++) = 0;
    *(scor_i++) = 0;
    *(slow_i++) = slow_energy;
    *(slow_i++) = slow_gradient;
    *(slow_i++) = 0;
    *(slow_i++) = 0;

  }

  // patch up data for i=0, in particular slow_
  fast_table[0] = fast_table[4] - fast_table[5] * r2_delta;
  fast_table[1] = fast_table[5];  // fast_gradient
  fast_table[2] = 0;
  fast_table[3] = 0;
  scor_table[0] = scor_table[4] - scor_table[5] * r2_delta;
  scor_table[1] = scor_table[5];  // scor_gradient
  scor_table[2] = 0;
  scor_table[3] = 0;
  slow_table[0] = slow_table[4] - slow_table[5] * r2_delta;
  slow_table[1] = slow_table[5];  // slow_gradient
  slow_table[2] = 0;
  slow_table[3] = 0;

  int j;
  for ( j=0; j<3; ++j ) {
    BigReal *t0 = 0;
    switch (j) {
      case 0: 
        t0 = fast_table;
      break;
      case 1: 
        t0 = scor_table;
      break;
      case 2: 
        t0 = slow_table;
      break;
    }
    BigReal *t;
    for ( i=0,t=t0; i<(n-1); ++i,t+=4 ) {
      BigReal x = r2_delta;
      BigReal v1 = t[0];
      BigReal g1 = t[1];
      BigReal v2 = t[4];
      BigReal g2 = t[5];
      // explicit formulas for v1 + g1 x + c x^2 + d x^3
      BigReal c = ( 3.0 * (v2 - v1) - x * (2.0 * g1 + g2) ) / ( x * x );
      BigReal d = ( -2.0 * (v2 - v1) + x * (g1 + g2) ) / ( x * x * x );
      // since v2 - v1 is imprecise, we refine c and d numerically
      // important because we need accurate forces (more than energies!)
      for ( int k=0; k < 2; ++k ) {
        BigReal dv = (v1 - v2) + ( ( d * x + c ) * x + g1 ) * x;
        BigReal dg = (g1 - g2) + ( 3.0 * d * x + 2.0 * c ) * x;
        c -= ( 3.0 * dv - x * dg ) / ( x * x );
        d -= ( -2.0 * dv + x * dg ) / ( x * x * x );
      }
      // store in the array;
      t[2] = c;  t[3] = d;
    }
    BigReal dvmax = 0;
    BigReal dgmax = 0;
    for ( i=0,t=t0; i<(n-1); ++i,t+=4 ) {
      BigReal x = r2_delta;
      BigReal dv = ( ( t[3] * x + t[2] ) * x + t[1] ) * x + t[0] - t[4];
      BigReal dg = ( 3.0 * t[3] * x + 2.0 * t[2] ) * x + t[1] - t[5];
      if ( fabs(dv) > dvmax ) dvmax = fabs(dv);
      if ( fabs(dg) > dgmax ) dgmax = fabs(dg);
      // if ( dv != 0.0 ) CkPrintf("TABLE %d ENERGY ERROR %g AT %g\n",j,dv,x*i);
      // if ( dg != 0.0 ) CkPrintf("TABLE %d FORCE ERROR %g AT %g\n",j,dg,x*i);
    }
    if ( ( ( dvmax != 0.0 ) || ( dgmax != 0.0 ) ) && ! CkMyPe() ) {
      iout << iINFO << "NONZERO IMPRECISION IN COULOMB TABLE: " <<
				dvmax << " " << dgmax << "\n" << endi;
    }
  }

#if 0
  char fname[100];
  sprintf(fname,"/tmp/namd.table.pe%d.dat",CkMyPe());
  FILE *f = fopen(fname,"w");
  for ( i=0; i<(n-1); ++i ) {
    BigReal r2 = r2_delta * i;
    BigReal *t;
    fprintf(f,"%g",r2);
    t = fast_table + 4*i;
    fprintf(f," %g %g %g %g", t[0], t[1], t[2], t[3]);
    t = scor_table + 4*i;
    fprintf(f," %g %g %g %g", t[0], t[1], t[2], t[3]);
    t = slow_table + 4*i;
    fprintf(f," %g %g %g %g", t[0], t[1], t[2], t[3]);
    fprintf(f,"\n");
  }
  fclose(f);
#endif

}

// clear all
// define interaction type (pair or self)
#define NBPAIR	1
#define NBSELF	2

#define NBTYPE NBPAIR
#include "ComputeNonbondedBase.h"
#define FULLELECT
#include "ComputeNonbondedBase.h"
#define SLOWONLY
#include "ComputeNonbondedBase.h"
#undef SLOWONLY
#undef FULLELECT
#undef  NBTYPE

#define NBTYPE NBSELF
#include "ComputeNonbondedBase.h"
#define FULLELECT
#include "ComputeNonbondedBase.h"
#define SLOWONLY
#include "ComputeNonbondedBase.h"
#undef SLOWONLY
#undef FULLELECT
#undef  NBTYPE

