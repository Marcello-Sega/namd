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

typedef unsigned short plint;

class Pairlists {
  enum {initsize = 10};
  plint *data;
  int curpos;
  int size;
  Pairlists(const Pairlists&) { ; }
  Pairlists& operator=(const Pairlists&) { return *this; }
public:
  Pairlists() : size(initsize) { data = new plint[initsize]; }
  ~Pairlists() { delete [] data; }
  plint *newlist(int max_size) {  // get a new list w/ room for max_size
    int reqnewsize = curpos + max_size + 1;
    int newsize = size;
    while ( newsize < reqnewsize ) { newsize += newsize >> 1; }
    if ( newsize > size ) {
      plint *newdata = new plint[newsize];
      CmiMemcpy(newdata,data,curpos*sizeof(plint));
      delete [] data;
      data = newdata;
      size = newsize;
    }
    return &data[curpos+1];
  }
  void newsize(int list_size) {  // set the size of the last list gotten
    data[curpos] = list_size;
    curpos += list_size + 1;
  }
  void reset() { curpos = 0; }  // go back to the beginning
  void nextlist(plint **list, int *list_size) {  // get next list and size
    *list = &data[curpos+1];
    curpos += ( *list_size = data[curpos] ) + 1;
  }
};

#define NBWORKARRAYSINIT(ARRAYS) \
  ComputeNonbondedWorkArrays* const computeNonbondedWorkArrays = ARRAYS;

#define NBWORKARRAY(TYPE,NAME,SIZE) \
  computeNonbondedWorkArrays->NAME.resize(SIZE); \
  TYPE * NAME = computeNonbondedWorkArrays->NAME.begin();

class ComputeNonbondedWorkArrays {
public:
  ResizeArray<plint> pairlisti;
  ResizeArray<BigReal> r2list;
  ResizeArray<plint> grouplist;
  ResizeArray<plint> fixglist;
  ResizeArray<plint> goodglist;
  ResizeArray<plint> pairlistx;
  ResizeArray<plint> pairlistm;
  ResizeArray<plint> pairlist;
  ResizeArray<plint> pairlist2;
  ResizeArray<short> vdwtype_array;
  ResizeArray<Force> f_0;
  ResizeArray<Force> fullf_0;
};

// function arguments
struct nonbonded {
  CompAtom* p[2];
#ifdef MEM_OPT_VERSION
  CompAtomExt *pExt[2];
#endif
  Force* ff[2];
  // for full electrostatics
  Force* fullf [2];

  int numAtoms[2];

  BigReal *reduction;
  BigReal *pressureProfileReduction;

  ComputeNonbondedWorkArrays *workArrays;

  Pairlists *pairlists;
  int savePairlists;
  int usePairlists;
  BigReal plcutoff;
  BigReal groupplcutoff;

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
  static void (*calcPairEnergy)(nonbonded *);
  static void (*calcSelf)(nonbonded *);
  static void (*calcSelfEnergy)(nonbonded *);

  static void (*calcFullPair)(nonbonded *);
  static void (*calcFullPairEnergy)(nonbonded *);
  static void (*calcFullSelf)(nonbonded *);
  static void (*calcFullSelfEnergy)(nonbonded *);

  static void (*calcMergePair)(nonbonded *);
  static void (*calcMergePairEnergy)(nonbonded *);
  static void (*calcMergeSelf)(nonbonded *);
  static void (*calcMergeSelfEnergy)(nonbonded *);

  static void (*calcSlowPair)(nonbonded *);
  static void (*calcSlowPairEnergy)(nonbonded *);
  static void (*calcSlowSelf)(nonbonded *);
  static void (*calcSlowSelfEnergy)(nonbonded *);

  enum { exclChecksumIndex, pairlistWarningIndex,
	 electEnergyIndex, fullElectEnergyIndex, vdwEnergyIndex,
//sd-db
	 electEnergyIndex_s, fullElectEnergyIndex_s, vdwEnergyIndex_s,
//sd-de
	 TENSOR(virialIndex), TENSOR(fullElectVirialIndex),
         VECTOR(pairVDWForceIndex), VECTOR(pairElectForceIndex),
	 reductionDataSize };
  static void submitReductionData(BigReal*,SubmitReduction*);
  static void submitPressureProfileData(BigReal*,SubmitReduction*);

  static Bool commOnly;
  static Bool fixedAtomsOn;
  static BigReal cutoff;
  static BigReal cutoff2;
  static BigReal dielectric_1;
  static const LJTable* ljTable;
  static const Molecule* mol;
  static BigReal r2_delta, r2_delta_1;
  static int r2_delta_exp;
  static BigReal *table_alloc;
  static BigReal *table_short;
  static BigReal *table_noshort;
  static BigReal *fast_table;
  static BigReal *scor_table;
  static BigReal *slow_table;
  static BigReal *corr_table;
  static BigReal *full_table;
  static BigReal *vdwa_table;
  static BigReal *vdwb_table;
  static BigReal *r2_table;
  static BigReal scaling;
  static BigReal scale14;
  static BigReal switchOn;
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
//sd-db
  static Bool fepOn;
  static BigReal lambda;
  static BigReal lambda2;
  static BigReal fepVdwShiftCoeff;
  static BigReal fepElecLambdaStart;
  static BigReal fepVdwLambdaEnd;
//sd-de
  static Bool lesOn;
  static int lesFactor;
  static BigReal lesScaling;

  static BigReal *lambda_table;

  static Bool pairInteractionOn;
  static Bool pairInteractionSelf;

  static Bool pressureProfileOn;
  static int pressureProfileSlabs;
  static int pressureProfileAtomTypes;
  static BigReal pressureProfileThickness;
  static BigReal pressureProfileMin;

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

  static void calc_error(nonbonded *);

  static void calc_pair(nonbonded *);
  static void calc_pair_energy(nonbonded *);
  static void calc_pair_fullelect(nonbonded *);
  static void calc_pair_energy_fullelect(nonbonded *);
  static void calc_pair_merge_fullelect(nonbonded *);
  static void calc_pair_energy_merge_fullelect(nonbonded *);
  static void calc_pair_slow_fullelect(nonbonded *);
  static void calc_pair_energy_slow_fullelect(nonbonded *);

  static void calc_self(nonbonded *);
  static void calc_self_energy(nonbonded *);
  static void calc_self_fullelect(nonbonded *);
  static void calc_self_energy_fullelect(nonbonded *);
  static void calc_self_merge_fullelect(nonbonded *);
  static void calc_self_energy_merge_fullelect(nonbonded *);
  static void calc_self_slow_fullelect(nonbonded *);
  static void calc_self_energy_slow_fullelect(nonbonded *);

//alchemical fep calcualtion
  static void calc_pair_energy_fep(nonbonded *);
  static void calc_pair_energy_fullelect_fep (nonbonded *);
  static void calc_pair_energy_merge_fullelect_fep (nonbonded *);
  static void calc_pair_energy_slow_fullelect_fep (nonbonded *);
  static void calc_self_energy_fep (nonbonded *);
  static void calc_self_energy_fullelect_fep (nonbonded *);
  static void calc_self_energy_merge_fullelect_fep (nonbonded *);
  static void calc_self_energy_slow_fullelect_fep (nonbonded *);

//locally enhanced sampling calcualtion
  static void calc_pair_les(nonbonded *);
  static void calc_pair_energy_les(nonbonded *);
  static void calc_pair_fullelect_les (nonbonded *);
  static void calc_pair_energy_fullelect_les (nonbonded *);
  static void calc_pair_merge_fullelect_les (nonbonded *);
  static void calc_pair_energy_merge_fullelect_les (nonbonded *);
  static void calc_pair_slow_fullelect_les (nonbonded *);
  static void calc_pair_energy_slow_fullelect_les (nonbonded *);
  static void calc_self_les (nonbonded *);
  static void calc_self_energy_les (nonbonded *);
  static void calc_self_fullelect_les (nonbonded *);
  static void calc_self_energy_fullelect_les (nonbonded *);
  static void calc_self_merge_fullelect_les (nonbonded *);
  static void calc_self_energy_merge_fullelect_les (nonbonded *);
  static void calc_self_slow_fullelect_les (nonbonded *);
  static void calc_self_energy_slow_fullelect_les (nonbonded *);

//pair interaction calcualtion
  static void calc_pair_energy_int(nonbonded *);
  static void calc_pair_energy_fullelect_int(nonbonded *);
  static void calc_pair_energy_merge_fullelect_int(nonbonded *);
  static void calc_self_energy_int (nonbonded *);
  static void calc_self_energy_fullelect_int(nonbonded *);
  static void calc_self_energy_merge_fullelect_int(nonbonded *);

//pressure profile calcualtion
  static void calc_pair_pprof(nonbonded *);
  static void calc_pair_energy_pprof(nonbonded *);
  static void calc_pair_fullelect_pprof(nonbonded *);
  static void calc_pair_energy_fullelect_pprof(nonbonded *);
  static void calc_pair_merge_fullelect_pprof(nonbonded *);
  static void calc_pair_energy_merge_fullelect_pprof(nonbonded *);
  static void calc_pair_slow_fullelect_pprof(nonbonded *);
  static void calc_pair_energy_slow_fullelect_pprof(nonbonded *);
  static void calc_self_pprof(nonbonded *);
  static void calc_self_energy_pprof(nonbonded *);
  static void calc_self_fullelect_pprof (nonbonded *);
  static void calc_self_energy_fullelect_pprof (nonbonded *);
  static void calc_self_merge_fullelect_pprof (nonbonded *);
  static void calc_self_energy_merge_fullelect_pprof (nonbonded *);
  static void calc_self_slow_fullelect_pprof (nonbonded *);
  static void calc_self_energy_slow_fullelect_pprof (nonbonded *);
};

#endif

