/**
***  Copyright (c) 1995, 1996, 1997, 1998, 1999, 2000 by
***  The Board of Trustees of the University of Illinois.
***  All rights reserved.
**/

/*
   Common operations for ComputeNonbonded classes
*/

#ifndef COMPUTENONBONDEDUTIL_H
#define COMPUTENONBONDEDUTIL_H

#include "NamdTypes.h"
#include "ReductionMgr.h"
class LJTable;
class Molecule;

// function arguments
struct nonbonded {
  CompAtom* p[2];
  Force* ff[2];
  // for full electrostatics
  Force* fullf [2];

  // used by pair and self
  int numAtoms[2];
  BigReal *reduction;

  // used by excl
  Position p_ij;
  int m14;

  // used by self
  int minPart;
  int maxPart;
  int numParts;

};

class ComputeNonbondedUtil {

public:

  ComputeNonbondedUtil() {}
  virtual ~ComputeNonbondedUtil() {}
  static void select(void);

  static void (*calcPair)(nonbonded *);
  static void (*calcSelf)(nonbonded *);
  static void (*calcExcl)(nonbonded *);

  static void (*calcFullPair)(nonbonded *);
  static void (*calcFullSelf)(nonbonded *);
  static void (*calcFullExcl)(nonbonded *);

  static void (*calcSlowPair)(nonbonded *);
  static void (*calcSlowSelf)(nonbonded *);
  static void (*calcSlowExcl)(nonbonded *);

  enum { exclChecksumIndex,
	 electEnergyIndex, fullElectEnergyIndex, vdwEnergyIndex,
	 TENSOR(virialIndex), TENSOR(fullElectVirialIndex),
	 reductionDataSize };
  static void submitReductionData(BigReal*,SubmitReduction*);

  static Bool commOnly;
  static Bool fixedAtomsOn;
  static Real cutoff;
  static BigReal cutoff2;
  static BigReal groupcutoff2;
  static BigReal dielectric_1;
  static const LJTable* ljTable;
  static const Molecule* mol;
  static BigReal scaling;
  static BigReal scale14;
  static Real switchOn;
  static BigReal switchOn_1;
  static BigReal switchOn2;
  static BigReal c0;
  static BigReal c1;
  static BigReal c3;
  static BigReal c5;
  static BigReal c6;
  static BigReal c7;
  static BigReal c8;
  // static BigReal d0;

  // for particle mesh Ewald
  static BigReal ewaldcof;
  static BigReal pi_ewaldcof;

  // for simplifying some common functions
  inline static BigReal square(	const BigReal &x,
 				const BigReal &y,
				const BigReal &z)
	{
	return(x*x+y*y+z*z);
	}

  // for splitting
  // No splitting
  static void calc_pair (nonbonded *);
  static void calc_pair_fullelect (nonbonded *);
  static void calc_pair_fullelect_pme (nonbonded *);
  static void calc_pair_slow_fullelect (nonbonded *);
  static void calc_pair_slow_fullelect_pme (nonbonded *);
  static void calc_self (nonbonded *);
  static void calc_self_fullelect (nonbonded *);
  static void calc_self_fullelect_pme (nonbonded *);
  static void calc_self_slow_fullelect (nonbonded *);
  static void calc_self_slow_fullelect_pme (nonbonded *);
  static void calc_excl (nonbonded *);
  static void calc_excl_fullelect (nonbonded *);
  static void calc_excl_fullelect_pme (nonbonded *);
  static void calc_excl_slow_fullelect (nonbonded *);
  static void calc_excl_slow_fullelect_pme (nonbonded *);
  inline static void shifting(BigReal &, BigReal &,
		       const BigReal &, const BigReal &,
		       const BigReal &, const BigReal &);

  // C1 Splitting
  static void calc_pair_c1 (nonbonded *);
  static void calc_pair_fullelect_c1 (nonbonded *);
  static void calc_pair_fullelect_pme_c1 (nonbonded *);
  static void calc_pair_slow_fullelect_c1 (nonbonded *);
  static void calc_pair_slow_fullelect_pme_c1 (nonbonded *);
  static void calc_self_c1 (nonbonded *);
  static void calc_self_fullelect_c1 (nonbonded *);
  static void calc_self_fullelect_pme_c1 (nonbonded *);
  static void calc_self_slow_fullelect_c1 (nonbonded *);
  static void calc_self_slow_fullelect_pme_c1 (nonbonded *);
  static void calc_excl_c1 (nonbonded *);
  static void calc_excl_fullelect_c1 (nonbonded *);
  static void calc_excl_fullelect_pme_c1 (nonbonded *);
  static void calc_excl_slow_fullelect_c1 (nonbonded *);
  static void calc_excl_slow_fullelect_pme_c1 (nonbonded *);
  inline static void c1splitting(BigReal &, BigReal &,
 			  const BigReal &, const BigReal &, const BigReal &);

  // XPLOR Splitting
  static void calc_pair_xplor (nonbonded *);
  static void calc_pair_fullelect_xplor (nonbonded *);
  static void calc_pair_fullelect_pme_xplor (nonbonded *);
  static void calc_pair_slow_fullelect_xplor (nonbonded *);
  static void calc_pair_slow_fullelect_pme_xplor (nonbonded *);
  static void calc_self_xplor (nonbonded *);
  static void calc_self_fullelect_xplor (nonbonded *);
  static void calc_self_fullelect_pme_xplor (nonbonded *);
  static void calc_self_slow_fullelect_xplor (nonbonded *);
  static void calc_self_slow_fullelect_pme_xplor (nonbonded *);
  static void calc_excl_xplor (nonbonded *);
  static void calc_excl_fullelect_xplor (nonbonded *);
  static void calc_excl_fullelect_pme_xplor (nonbonded *);
  static void calc_excl_slow_fullelect_xplor (nonbonded *);
  static void calc_excl_slow_fullelect_pme_xplor (nonbonded *);
  inline static void xplorsplitting(BigReal &, BigReal &,
			     const BigReal &,const BigReal &);

};

#endif

