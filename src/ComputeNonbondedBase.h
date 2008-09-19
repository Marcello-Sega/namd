/**
***  Copyright (c) 1995, 1996, 1997, 1998, 1999, 2000 by
***  The Board of Trustees of the University of Illinois.
***  All rights reserved.
**/

// Several special cases are defined:
//   NBTYPE: exclusion method (NBPAIR, NBSELF -- mutually exclusive)
//   FULLELECT full electrostatics calculation?

#ifdef ARCH_POWERPC
#include <builtins.h>
#endif

#if defined(NAMD_SSE) && defined(__INTEL_COMPILER) && defined(__SSE2__)
#include <emmintrin.h>  // We're using SSE2 intrinsics
#endif

#ifdef DEFINITION // (
  #include "LJTable.h"
  #include "Molecule.h"
  #include "ComputeNonbondedUtil.h"
#endif // )

// determining class name
#undef NAME
#undef CLASS
#undef CLASSNAME
#define NAME CLASSNAME(calc)

#undef PAIR
#if NBTYPE == NBPAIR
  #define PAIR(X) X
  #define CLASS ComputeNonbondedPair
  #define CLASSNAME(X) ENERGYNAME( X ## _pair )
#else
  #define PAIR(X)
#endif

#undef SELF
#if NBTYPE == NBSELF
  #define SELF(X) X
  #define CLASS ComputeNonbondedSelf
  #define CLASSNAME(X) ENERGYNAME( X ## _self )
#else
  #define SELF(X)
#endif

#undef ENERGYNAME
#undef ENERGY
#undef NOENERGY
#ifdef CALCENERGY
  #define ENERGY(X) X
  #define NOENERGY(X)
  #define ENERGYNAME(X) SLOWONLYNAME( X ## _energy )
#else
  #define ENERGY(X)
  #define NOENERGY(X) X
  #define ENERGYNAME(X) SLOWONLYNAME( X )
#endif

#undef SLOWONLYNAME
#undef FAST
#ifdef SLOWONLY
  #define FAST(X)
  #define SLOWONLYNAME(X) MERGEELECTNAME( X ## _slow )
#else
  #define FAST(X) X
  #define SLOWONLYNAME(X) MERGEELECTNAME( X )
#endif

#undef MERGEELECTNAME
#undef SHORT
#undef NOSHORT
#ifdef MERGEELECT
  #define SHORT(X)
  #define NOSHORT(X) X
  #define MERGEELECTNAME(X) FULLELECTNAME( X ## _merge )
#else
  #define SHORT(X) X
  #define NOSHORT(X)
  #define MERGEELECTNAME(X) FULLELECTNAME( X )
#endif

#undef FULLELECTNAME
#undef FULL
#undef NOFULL
#ifdef FULLELECT
  #define FULLELECTNAME(X) FEPNAME( X ## _fullelect )
  #define FULL(X) X
  #define NOFULL(X)
#else
  #define FULLELECTNAME(X) FEPNAME( X )
  #define FULL(X)
  #define NOFULL(X) X
#endif

// Here are what these macros stand for:
// FEP/NOT_FEP: FEP free energy perturbation is active/inactive 
//      (does NOT use LAM)
// LES: locally-enhanced sampling is active
// LAM: scale the potential by a factor lambda (except FEP)
// INT: measure interaction energies
// PPROF: pressure profiling

#undef FEPNAME
#undef FEP
#undef LES
#undef INT
#undef PPROF
#undef LAM
#undef ALCH
#undef TI
#define FEPNAME(X) LAST( X )
#define FEP(X)
#define ALCHPAIR(X)
#define NOT_ALCHPAIR(X) X
#define LES(X)
#define INT(X)
#define PPROF(X)
#define LAM(X)
#define ALCH(X)
#define TI(X)
#ifdef FEPFLAG
  #undef FEPNAME
  #undef FEP
  #undef ALCH
  #define FEPNAME(X) LAST( X ## _fep )
  #define FEP(X) X
  #define ALCH(X) X
#endif
#ifdef TIFLAG
  #undef FEPNAME
  #undef TI
  #undef ALCH
  #define FEPNAME(X) LAST( X ## _ti )
  #define TI(X) X
  #define ALCH(X) X
#endif
#ifdef LESFLAG
  #undef FEPNAME
  #undef LES
  #undef LAM
  #define FEPNAME(X) LAST( X ## _les )
  #define LES(X) X
  #define LAM(X) X
#endif
#ifdef INTFLAG
  #undef FEPNAME
  #undef INT
  #define FEPNAME(X) LAST( X ## _int )
  #define INT(X) X
#endif
#ifdef PPROFFLAG
  #undef FEPNAME
  #undef INT
  #undef PPROF
  #define FEPNAME(X) LAST( X ## _pprof )
  #define INT(X) X
  #define PPROF(X) X
#endif

#define LAST(X) X

// see if things are really messed up
SELF( PAIR( foo bar ) )
LES( FEP( foo bar ) )
LES( INT( foo bar ) )
FEP( INT( foo bar ) )
LAM( INT( foo bar ) )
FEP( NOENERGY( foo bar ) )
ENERGY( NOENERGY( foo bar ) )


// ************************************************************
// function header
void ComputeNonbondedUtil :: NAME
  ( nonbonded *params )

// function body
{
  // int NAME;  // to track errors via unused variable warnings

  if ( ComputeNonbondedUtil::commOnly ) return;

  // speedup variables
  BigReal *reduction = params->reduction;

  PPROF(
  BigReal *pressureProfileReduction = params->pressureProfileReduction;
  const BigReal invThickness = 1.0 / pressureProfileThickness;
  )

  Pairlists &pairlists = *(params->pairlists);
  int savePairlists = params->savePairlists;
  int usePairlists = params->usePairlists;
  pairlists.reset();
  // PAIR(iout << "--------\n" << endi;)

  // local variables
  int exclChecksum = 0;
  FAST
  (
  ENERGY( BigReal vdwEnergy = 0; )
  SHORT
  (
  ENERGY( BigReal electEnergy = 0; )
  )

  FEP
  (
  ENERGY( BigReal vdwEnergy_s = 0; )
  SHORT
  (
  ENERGY( BigReal electEnergy_s = 0; )
  )
  )
  
  SHORT
  (
  BigReal virial_xx = 0;
  BigReal virial_xy = 0;
  BigReal virial_xz = 0;
  BigReal virial_yy = 0;
  BigReal virial_yz = 0;
  BigReal virial_zz = 0;
  )
  )
  FULL
  (
  ENERGY( BigReal fullElectEnergy = 0; )
  FEP
  (
  ENERGY( BigReal fullElectEnergy_s = 0; )
  )
  BigReal fullElectVirial_xx = 0;
  BigReal fullElectVirial_xy = 0;
  BigReal fullElectVirial_xz = 0;
  BigReal fullElectVirial_yy = 0;
  BigReal fullElectVirial_yz = 0;
  BigReal fullElectVirial_zz = 0;
  )

  // Bringing stuff into local namespace for speed.

  const BigReal offset_x = params->offset.x;
  const BigReal offset_y = params->offset.y;
  const BigReal offset_z = params->offset.z;

  register const BigReal plcutoff2 = \
 			params->plcutoff * params->plcutoff;
  register const BigReal groupplcutoff2 = \
	 		params->groupplcutoff * params->groupplcutoff;
  const BigReal dielectric_1 = ComputeNonbondedUtil:: dielectric_1;
  const LJTable* const ljTable = ComputeNonbondedUtil:: ljTable;
  LJTable::TableEntry ljNull;  ljNull.A = 0; ljNull.B = 0;
  const LJTable::TableEntry* const lj_null_pars = &ljNull;
  const Molecule* const mol = ComputeNonbondedUtil:: mol;
  SHORT
  (
  const BigReal* const table_four = ComputeNonbondedUtil:: table_short;
  )
  FULL
  (
  SHORT
  (
  const BigReal* const slow_table = ComputeNonbondedUtil:: slow_table;
  )
  NOSHORT
  (
//#if 1 ALCH(-1)
  const BigReal* const table_four = ComputeNonbondedUtil:: table_noshort;
//#else  // have to switch this for ALCH
//  BigReal* table_four = ComputeNonbondedUtil:: table_noshort;
//#endif
  )
  )
  const BigReal scaling = ComputeNonbondedUtil:: scaling;
  const BigReal modf_mod = 1.0 - scale14;
  FAST
  (
  const BigReal switchOn2 = ComputeNonbondedUtil:: switchOn2;
  const BigReal c1 = ComputeNonbondedUtil:: c1;
  const BigReal c3 = ComputeNonbondedUtil:: c3;
  )
  const BigReal r2_delta = ComputeNonbondedUtil:: r2_delta;
  const int r2_delta_exp = ComputeNonbondedUtil:: r2_delta_exp;
  // const int r2_delta_expc = 64 * (r2_delta_exp - 127);
  const int r2_delta_expc = 64 * (r2_delta_exp - 1023);

  ALCH(
    const BigReal switchdist2 = ComputeNonbondedUtil::switchOn2;
    const BigReal cutoff2 = ComputeNonbondedUtil::cutoff2;
    const BigReal switchfactor = 1./((cutoff2 - switchdist2)*(cutoff2 - switchdist2)*(cutoff2 - switchdist2));
    const BigReal fepElecLambdaStart = ComputeNonbondedUtil::fepElecLambdaStart;
    const BigReal fepVdwLambdaEnd = ComputeNonbondedUtil::fepVdwLambdaEnd;
    const BigReal fepVdwShiftCoeff = ComputeNonbondedUtil::fepVdwShiftCoeff;

    /*lambda values 'up' are for atoms scaled up with lambda (partition 1)*/
    BigReal lambdaUp = ComputeNonbondedUtil::lambda;
    BigReal elecLambdaUp =  (lambdaUp <= fepElecLambdaStart)? 0. : \
              (lambdaUp - fepElecLambdaStart) / (1. - fepElecLambdaStart);
    BigReal vdwLambdaUp = 
        (lambdaUp >= fepVdwLambdaEnd)? 1. : lambdaUp / fepVdwLambdaEnd; 
    BigReal vdwShiftUp = ComputeNonbondedUtil::fepVdwShiftCoeff * (1-vdwLambdaUp);
    FEP(BigReal lambda2Up = ComputeNonbondedUtil::lambda2;)
    FEP(BigReal elecLambda2Up = (lambda2Up <= fepElecLambdaStart)? 0. : \
              (lambda2Up - fepElecLambdaStart) / (1. - fepElecLambdaStart);)
    FEP(BigReal vdwLambda2Up = 
        (lambda2Up >= fepVdwLambdaEnd)? 1. : lambda2Up / fepVdwLambdaEnd;) 
    FEP(BigReal vdwShift2Up = ComputeNonbondedUtil::fepVdwShiftCoeff * (1-vdwLambda2Up);)

        
    /*lambda values 'down' are for atoms scaled down with lambda (partition 2)*/
    BigReal lambdaDown = 1 - ComputeNonbondedUtil::lambda;
    BigReal elecLambdaDown =  (lambdaDown <= fepElecLambdaStart)? 0. : \
              (lambdaDown - fepElecLambdaStart) / (1. - fepElecLambdaStart);
    BigReal vdwLambdaDown = 
        (lambdaDown >= fepVdwLambdaEnd)? 1. : lambdaDown / fepVdwLambdaEnd; 
    BigReal vdwShiftDown = ComputeNonbondedUtil::fepVdwShiftCoeff * (1-vdwLambdaDown);
    FEP(BigReal lambda2Down = 1 - ComputeNonbondedUtil::lambda2;)
    FEP(BigReal elecLambda2Down = (lambda2Down <= fepElecLambdaStart)? 0. : \
              (lambda2Down - fepElecLambdaStart) / (1. - fepElecLambdaStart); )
    FEP(BigReal vdwLambda2Down = 
        (lambda2Down >= fepVdwLambdaEnd)? 1. : lambda2Down / fepVdwLambdaEnd; )
    FEP(BigReal vdwShift2Down = ComputeNonbondedUtil::fepVdwShiftCoeff * (1-vdwLambda2Down);)

  
  // Thermodynamic Integration Notes (F. Buelens): 
  // Separation of atom pairs into different pairlist according to lambda
  // coupling, for alchemical free energy calculations. Normal pairlists
  // (pairlist[nxm]_save) are re-used for non-lambda-coupled pairs; new ones
  // (pairlist[nxm][12]_save are created for atoms switched up or down with
  // lambda respectively.
  // This makes TI coding far easier and more readable, since it's necessary 
  // to store dU/dlambda in different variables depending on whether atoms are
  // being switched up or down. Further, it allows Jordi's explicit coding of 
  // the separation-shifted vdW potential to stay in place without a big 
  // performance hit, since it's no longer necessary to calculate vdW potentials
  // explicity for the bulk (non-alchemical) interactions - so that part of the 
  // free energy code stays readable and easy to modify. Finally there should
  // be some performance gain over the old FEP implementation as we no
  // longer have to check the partitions of each atom pair and switch
  // parameters accordingly for every single NonbondedBase2.h loop - this is 
  // done at the pairlist level
  
  int pswitchTable[3*3];
  // determines which pairlist alchemical pairings get sent to
  // 0: non-alchemical pairs, partition 0 <-> partition 0
  // 1: atoms scaled up as lambda increases, p0<->p1
  // 2: atoms scaled down as lambda increases, p0<->p2
  // all p1<->p2 interactions to be dropped (99)
  // in general, 'up' refers to 1, 'down' refers to 2
  for (int ip=0; ip<3; ++ip) {
    for (int jp=0; jp<3; ++jp ) {
      pswitchTable[ip+3*jp] = 0;
      if ((ip==1 && jp==0) || (ip==0 && jp==1)) pswitchTable[ip+3*jp] = 1;
      if ((ip==2 && jp==0) || (ip==0 && jp==2)) pswitchTable[ip+3*jp] = 2;
      if (ip + jp == 3) pswitchTable[ip+3*jp] = 99; // no interaction

      if (! ComputeNonbondedUtil::decouple) {
        // no decoupling: interactions within a partition are treated like
        // normal alchemical pairs
        if (ip == 1 && jp == 1) pswitchTable[ip+3*jp] = 1;
        if (ip == 2 && jp == 2) pswitchTable[ip+3*jp] = 2;
      }
      if (ComputeNonbondedUtil::decouple) {
        // decoupling: PME calculates extra grids so that while PME 
        // interaction with the full system is switched off, a new PME grid
        // containing only alchemical atoms is switched on. Full interactions 
        // between alchemical atoms are maintained; potentials within one 
        // partition need not be scaled here.
        if (ip == 1 && jp == 1) pswitchTable[ip+3*jp] = 0;
        if (ip == 2 && jp == 2) pswitchTable[ip+3*jp] = 0;
      }
    }
  }

  // dU/dlambda variables for thermodynamic integration
  TI(
      BigReal vdwEnergy_ti_1 = 0;
      BigReal vdwEnergy_ti_2 = 0;
      SHORT(BigReal electEnergy_ti_1 = 0;
      BigReal electEnergy_ti_2 = 0;)
      FULL(BigReal fullElectEnergy_ti_1 = 0; 
      BigReal fullElectEnergy_ti_2 = 0;) 
   )
  )
        
        
  const int i_upper = params->numAtoms[0];
  register const int j_upper = params->numAtoms[1];
  register int j;
  register int i;
  const CompAtom *p_0 = params->p[0];
  const CompAtom *p_1 = params->p[1];
#ifdef MEM_OPT_VERSION
  const CompAtomExt *pExt_0 = params->pExt[0];
  const CompAtomExt *pExt_1 = params->pExt[1];
#endif

  char * excl_flags_buff = 0;
  const int32 * full_excl = 0;
  const int32 * mod_excl = 0;

  plint *pairlistn_save;  int npairn;
  plint *pairlistx_save;  int npairx;
  plint *pairlistm_save;  int npairm;

  ALCH(
  plint *pairlistnA1_save;  int npairnA1;
  plint *pairlistxA1_save;  int npairxA1;
  plint *pairlistmA1_save;  int npairmA1;
  plint *pairlistnA2_save;  int npairnA2;
  plint *pairlistxA2_save;  int npairxA2;
  plint *pairlistmA2_save;  int npairmA2;
  )

  NBWORKARRAYSINIT(params->workArrays);

  int arraysize = j_upper+5;

  NBWORKARRAY(plint,pairlisti,arraysize)
  NBWORKARRAY(BigReal,r2list,arraysize)

  union { double f; int32 i[2]; } byte_order_test;
  byte_order_test.f = 1.0;  // should occupy high-order bits only
  int32 *r2iilist = (int32*)r2list + ( byte_order_test.i[0] ? 0 : 1 );

  if ( ! ( savePairlists || ! usePairlists ) ) arraysize = 0;


  // DMK - Atom Sort
  // 
  // Basic Idea: Determine the line between the center of masses of
  //   the two patches.  Project and then sort the lists of atoms
  //   along this line.  Then, as the pairlists are being generated
  //   for the atoms in the first atom list, use the sorted
  //   list to only add atoms from the second list that are between
  //   +/- ~cutoff from the atoms position on the line.
  #if NAMD_ComputeNonbonded_SortAtoms != 0 && ( 0 PAIR( + 1 ) )

  // NOTE: For the first atom list, just the sort values themselves need to be
  // calculated (a BigReal vs. SortEntry data structure).  However, a second
  // SortEntry array is needed to perform the merge sort on the second list of
  // atoms.  Use the atomSort_[0/1]_sortValues__ arrays to sort the second
  // list of atoms and then use the left over array to make a version of the
  // list that only includes non-fixed groups (fixg).

  NBWORKARRAY(SortEntry, atomSort_0_sortValues__, arraysize);
  NBWORKARRAY(SortEntry, atomSort_1_sortValues__, arraysize);
  NBWORKARRAY(BigReal, p_0_sortValues, arraysize);

  register SortEntry* p_1_sortValues = atomSort_0_sortValues__;
  register SortEntry* p_1_sortValues_fixg = atomSort_1_sortValues__;

  int p_0_sortValues_len = 0;
  int p_1_sortValues_len = 0;
  int p_1_sortValues_fixg_len = 0;

  BigReal atomSort_windowRadius = 0.0;

  if (savePairlists || !usePairlists) {

    /// Determine the two points to that will define the line ///

    // Center of mass for the first patch's atoms
    register BigReal p_0_avgX = 0.0;
    register BigReal p_0_avgY = 0.0;
    register BigReal p_0_avgZ = 0.0;
    {
      register int atomCount = 0;
      register int i = 0;

      register BigReal p_x = p_0->position.x;
      register BigReal p_y = p_0->position.y;
      register BigReal p_z = p_0->position.z;
      register int hgs = p_0->hydrogenGroupSize;

      while (i < i_upper) {

	i += hgs;
        register const CompAtom* p_i = p_0 + i;
        hgs = p_i->hydrogenGroupSize;

        p_0_avgX += p_x;
        p_0_avgY += p_y;
        p_0_avgZ += p_z;
        p_x = p_i->position.x;
        p_y = p_i->position.y;
        p_z = p_i->position.z;

        atomCount++;
      }

      p_0_avgX = p_0_avgX / ((double)atomCount);
      p_0_avgY = p_0_avgY / ((double)atomCount);
      p_0_avgZ = p_0_avgZ / ((double)atomCount);
    }

    // Center of mass for the second patch's atoms
    register BigReal p_1_avgX = 0.0;
    register BigReal p_1_avgY = 0.0;
    register BigReal p_1_avgZ = 0.0;
    {
      register int atomCount = 0;
      register int i = 0;

      register BigReal p_x = p_1->position.x;
      register BigReal p_y = p_1->position.y;
      register BigReal p_z = p_1->position.z;
      register int hgs = p_1->hydrogenGroupSize;

      while (i < j_upper) {

	i += hgs;
        register const CompAtom* p_i = p_1 + i;
        hgs = p_i->hydrogenGroupSize;

        p_1_avgX += p_x;
        p_1_avgY += p_y;
        p_1_avgZ += p_z;
        p_x = p_i->position.x;
        p_y = p_i->position.y;
        p_z = p_i->position.z;

        atomCount++;
      }

      p_1_avgX = p_1_avgX / ((double)atomCount);
      p_1_avgY = p_1_avgY / ((double)atomCount);
      p_1_avgZ = p_1_avgZ / ((double)atomCount);
    }

    // Need to move the points away from each other (so all of the
    // projected positions are between the two points).
    //   P0 = ( p_0_avgX, p_0_avgY, p_0_avgZ )
    //   P1 = ( p_1_avgX, p_1_avgY, p_1_avgZ )
    // P1 - P0 = V... so, V points from P0 towards P1
    //   Scale V so that |V| > corner to opposite corner distance of Patch
    //   NOTE: For now, scale to cutoff * 3
    register BigReal V_length;   // = |V|    (length of V)
    register BigReal V_length2;  // = |V|^2  (length of V squared)
    {
      register BigReal V_X = p_1_avgX - p_0_avgX;
      register BigReal V_Y = p_1_avgY - p_0_avgY;
      register BigReal V_Z = p_1_avgZ - p_0_avgZ;
      V_length2 = (V_X * V_X) + (V_Y * V_Y) + (V_Z * V_Z);
      V_length = sqrt(V_length2);
      register BigReal cutoff = ComputeNonbondedUtil::cutoff;
      register BigReal scaleFactor = (3.0 * cutoff) / V_length;
      V_X *= scaleFactor;
      V_Y *= scaleFactor;
      V_Z *= scaleFactor;

      // Move the points away from each other (P0 -= V, P1 += V)
      p_0_avgX -= V_X;
      p_0_avgY -= V_Y;
      p_0_avgZ -= V_Z;
      p_1_avgX += V_X;
      p_1_avgY += V_Y;
      p_1_avgZ += V_Z;

      // Recalculate V_length2 and V_length now that the points have moved
      V_length += 2.0 * (scaleFactor * V_length);
      V_length2 = V_length * V_length;
    }

    /// Create the list of sort values for each list of atoms (sort the second list) ///

    // For every atom, determine the sort value (i.e., the distance
    //   along 'the line, L, between P0 and P1' from P0 to the projection of
    //   the atom on L).
    // Let PA = position of atom in the patch space
    //   a = | P0 - PA |
    //   b = | P1 - PA |
    //   c = | P0 - P1 | = V_length
    // Case:  Triangle formed by P0, P1, and PA (position of atom).
    //   Projection of PA on the line L  between P0 and P1 (call this point
    //   PAp, note that line between P0 and P1 is perpendicular to the
    //   line between PA and PAp AND PAp is between P0 and P1 on L).
    // then ... 'distance between P0 and PAp' = (a^2 - b^2 + c^2) / (2 * c)
    // NOTE:  The '1 / (2 * c)' portion of this is constant since P0 and P1
    //   are now fixed.  Place this constant in a register to avoid repeating
    //   the divide calculation.
    register const BigReal V_length_multiplier = (1.0 / (2.0 * V_length));

    atomSort_windowRadius = params->groupplcutoff;

    // Atom list 1
    {
      register int i = 0;
      register unsigned int ngia = p_1->nonbondedGroupIsAtom;
      register unsigned int hgs = p_1->hydrogenGroupSize;
      register BigReal p_x = p_1->position.x;
      register BigReal p_y = p_1->position.y;
      register BigReal p_z = p_1->position.z;
      register int index = 0;

      // Try to take advantage of multiply-add instructions by pulling out the
      //   'V_length2 * V_length_multiplier' term in the sortVal calculation.
      //   So...  '(a2 - b2 + V_length2) * V_length_multiplier' becomes...
      //          '((a2 - b2) * V_length_multiplier) + (V_length2 * V_length_multiplier)'
      register const BigReal sortVal_addAmmount = V_length2 * V_length_multiplier;

      while (i < j_upper) {

	// Advance i... NOTE: ngia is either 0 or 1 so if ngia is set '1' is
	//   added to 'i'.  Otherwise, the value of 'hgs' is added to 'i'.
	i += (ngia) + ((1 - ngia) * hgs);

	// Update p_i_next to point to the atom for the next iteration and begin
	//   loading the 'ngia' and 'hgs' values for that atom.
        register const CompAtom* p_i_next = p_1 + i;
        ngia = p_i_next->nonbondedGroupIsAtom;
        hgs = p_i_next->hydrogenGroupSize;

        // Use the position values then start loading the next atom's position values.
        register BigReal a_x_diff = (p_0_avgX - p_x);
        register BigReal a_y_diff = (p_0_avgY - p_y);
        register BigReal a_z_diff = (p_0_avgZ - p_z);
        register BigReal b_x_diff = (p_1_avgX - p_x);
        register BigReal b_y_diff = (p_1_avgY - p_y);
        register BigReal b_z_diff = (p_1_avgZ - p_z);
        p_x = p_i_next->position.x;
        p_y = p_i_next->position.y;
        p_z = p_i_next->position.z;

        // Calculate and store the sort value for this atom (NOTE: c = V_length and c^2 = V_length2)
        register BigReal a2 = (a_x_diff * a_x_diff) + (a_y_diff * a_y_diff) + (a_z_diff * a_z_diff);
        register BigReal b2 = (b_x_diff * b_x_diff) + (b_y_diff * b_y_diff) + (b_z_diff * b_z_diff);
        register BigReal sortVal = (a2 - b2) * V_length_multiplier + sortVal_addAmmount;
	register SortEntry* p_1_sortValStorePtr = p_1_sortValues + p_1_sortValues_len;
        p_1_sortValStorePtr->index = index;
        p_1_sortValStorePtr->sortValue = sortVal;
        p_1_sortValues_len++;
        index = i;
      }
    }

    // NOTE: This list and another version of it with only non-fixed
    //   atoms will be used in place of grouplist and fixglist.
    {
      #if 0   // Selection Sort

        for (int i = 0; i < p_1_sortValues_len; i++) {

          // Search through the remaining elements, finding the lowest
          //   value, and then swap it with the first remaining element.
          //   Start by assuming the first element is the smallest.
          register int smallestIndex = i;
          register BigReal smallestValue = p_1_sortValues[i].sortValue;
          for (int j = i + 1; j < p_1_sortValues_len; j++) {
            register BigReal currentValue = p_1_sortValues[j].sortValue;
            if (currentValue < smallestValue) {
              smallestIndex = j;
              smallestValue = currentValue;
	    }
          }

          // Swap the first remaining element with the smallest element
          if (smallestIndex != i) {
            register SortEntry* entryA = p_1_sortValues + i;
            register SortEntry* entryB = p_1_sortValues + smallestIndex;
            register unsigned int tmpIndex = entryA->index;
            register BigReal tmpSortValue = entryA->sortValue;
            entryA->index = entryB->index;
            entryA->sortValue = entryB->sortValue;
            entryB->index = tmpIndex;
            entryB->sortValue = tmpSortValue;
	  }
        }

      #elif 0   // Bubble Sort

        register int keepSorting = 0;
        do {

          // Reset the keepSorting flag (assume no swaps will occur)
          keepSorting = 0;

          // Loop through the pairs and swap if needed
          register SortEntry* sortEntry1 = p_1_sortValues;
          for (int i = 1; i < p_1_sortValues_len; i++) {
            register SortEntry* sortEntry0 = sortEntry1;
            sortEntry1 = p_1_sortValues + i;
            register BigReal sortEntry0_sortValue = sortEntry0->sortValue;
            register BigReal sortEntry1_sortValue = sortEntry1->sortValue;
            if (sortEntry0_sortValue > sortEntry1_sortValue) {
              register int sortEntry0_index = sortEntry0->index;
              register int sortEntry1_index = sortEntry1->index;
              sortEntry0->index = sortEntry1_index;
              sortEntry0->sortValue = sortEntry1_sortValue;
              sortEntry1->index = sortEntry0_index;
              sortEntry1->sortValue = sortEntry0_sortValue;
              keepSorting = 1;
	    }
	  }

        } while (keepSorting != 0);  // Loop again if at least one set of
			             //   elements was swapped.
      #else   // Merge Sort

        #if NAMD_ComputeNonbonded_SortAtoms_LessBranches == 0

        register SortEntry* srcArray = p_1_sortValues;
        register SortEntry* dstArray = p_1_sortValues_fixg; //tmpArray;

        // Start with each element being a separate list.  Start
        //   merging the "lists" into larger lists.
        register int subListSize = 1;
        while (subListSize < p_1_sortValues_len) {

          // NOTE: This iteration consumes sublists of length
          //   subListSize and produces sublists of length
          //   (2*subListSize).  So, keep looping while the length of a
          //   single sorted sublist is not the size of the entire array.

          // Iterate through the lists, merging each consecutive pair of lists.
          register int firstListOffset = 0;
          while (firstListOffset < p_1_sortValues_len) {

            /// Setup pointers and counts for sublists in the pair. ///
            register int numElements = min(2 * subListSize, p_1_sortValues_len - firstListOffset);
            register int list0len;
            register int list1len;
            if (numElements > subListSize) {
              list0len = subListSize;                // First list full
              list1len = numElements - subListSize;  // 1+ elements in second list
	    } else {
              list0len = numElements;                // 1+ elements in first list
              list1len = 0;                          // Zero elements in second list
	    }

            register SortEntry* list0ptr = srcArray + firstListOffset;
            register SortEntry* list1ptr = list0ptr + subListSize;
            register SortEntry* dstptr = dstArray + firstListOffset;

            /// Merge the sublists ///

            // While there are elements in both lists, pick from one
            while (list0len > 0 && list1len > 0) {

              register BigReal sortValue0 = list0ptr->sortValue;
              register BigReal sortValue1 = list1ptr->sortValue;

              if (sortValue0 < sortValue1) {  // choose first list (list0)

                // Copy the value from srcArray to dstArray
                register int index0 = list0ptr->index;
                dstptr->sortValue = sortValue0;
                dstptr->index = index0;

                // Move the pointers forward for the sublists
                dstptr++;
                list0ptr++;
                list0len--;

              } else {                        // choose second list (list1)

                // Copy the value from srcArray to dstArray
                register int index1 = list1ptr->index;
                dstptr->sortValue = sortValue1;
                dstptr->index = index1;

                // Move the pointers forward for the sublists
                dstptr++;
                list1ptr++;
                list1len--;
              }

	    } // end while (list0len > 0 && list1len > 0)

            // NOTE: Either list0len or list1len is zero at this point
            //   so only one of the following loops should execute.

            // Drain remaining elements from the first list (list0)
            while (list0len > 0) {

              // Copy the value from srcArray to dstArray
              register BigReal sortValue0 = list0ptr->sortValue;
              register int index0 = list0ptr->index;
              dstptr->sortValue = sortValue0;
              dstptr->index = index0;

              // Move the pointers forward for the sublists
              dstptr++;
              list0ptr++;
              list0len--;

	    } // end while (list0len > 0)

            // Drain remaining elements from the first list (list1)
            while (list1len > 0) {

              // Copy the value from srcArray to dstArray
              register BigReal sortValue1 = list1ptr->sortValue;
              register int index1 = list1ptr->index;
              dstptr->sortValue = sortValue1;
              dstptr->index = index1;

              // Move the pointers forward for the sublists
              dstptr++;
              list1ptr++;
              list1len--;

	    } // end while (list1len > 0)

            // Move forward to the next pair of sub-lists
            firstListOffset += (2 * subListSize);

	  } // end while (firstListOffset < p_1_sortValues_len) {

          // Swap the dstArray and srcArray pointers
          register SortEntry* tmpPtr = dstArray;
          dstArray = srcArray;
          srcArray = tmpPtr;

          // Double the subListSize
          subListSize <<= 1;

	}  // end while (subListSize < p_1_sortValues_len)

        // Set the sort values pointers (NOTE: srcArray and dstArray are
        //   swapped at the end of each iteration of the merge sort outer-loop).
        p_1_sortValues_fixg = dstArray;
        p_1_sortValues = srcArray;

        #else

        // NOTE: This macro "returns" either val0 (if test == 0) or val1 (if
        // test == 1).  It expects test to be either 0 or 1 (no other values).
        #define TERNARY_ASSIGN(test, val0, val1)   ((test * val0) + ((1 - test) * val1))

        register SortEntry* srcArray = p_1_sortValues;
        register SortEntry* dstArray = p_1_sortValues_fixg; //tmpArray;

        // Start with each element being a separate list.  Start
        //   merging the "lists" into larger lists.
        register int subListSize = 1;
        while (subListSize < p_1_sortValues_len) {

          // NOTE: This iteration consumes sublists of length
          //   subListSize and produces sublists of length
          //   (2*subListSize).  So, keep looping while the length of a
          //   single sorted sublist is not the size of the entire array.

          // Iterate through the lists, merging each consecutive pair of lists.
          register int firstListOffset = 0;
          while (firstListOffset < p_1_sortValues_len) {

            /// Setup pointers and counts for sublists in the pair. ///

            // Calculate the number of elements for both sublists...
            //   min(2 * subListSize, p_1_sortValues_len - firstListOffset);
            register int numElements;
	    {
              register int numElements_val0 = 2 * subListSize;
              register int numElements_val1 = p_1_sortValues_len - firstListOffset;
              register bool numElements_test = (numElements_val0 < numElements_val1);
              numElements = TERNARY_ASSIGN(numElements_test, numElements_val0, numElements_val1);
	    }

            // Setup the pointers for the source and destination arrays
            register SortEntry* dstptr = dstArray + firstListOffset;    // destination array pointer
            register SortEntry* list0ptr = srcArray + firstListOffset;  // source list 0 pointer
            register SortEntry* list1ptr = list0ptr + subListSize;      // source list 1 pointer
            register SortEntry* list0ptr_end;  // pointer to end of source list0's elements (element after last)
            register SortEntry* list1ptr_end;  // pointer to end of source list1's elements (element after last)
            {
              register bool lenTest = (numElements > subListSize);
              register int list0len_val0 = subListSize;
              register int list1len_val0 = numElements - subListSize;
              register int list0len_val1 = numElements;  // NOTE: list1len_val1 = 0
              register int list0len = TERNARY_ASSIGN(lenTest, list0len_val0, list0len_val1);
              register int list1len = TERNARY_ASSIGN(lenTest, list1len_val0, 0);
              list0ptr_end = list0ptr + list0len;
              list1ptr_end = list1ptr + list1len;
	    }

            // The firstListOffset variable won't be used again until the next
            //   iteration, so go ahead and update it now...
            //   Move forward to the next pair of sub-lists
            firstListOffset += (2 * subListSize);

            /// Merge the sublists ///

            // Pre-load values from both source arrays
            register BigReal sortValue0 = list0ptr->sortValue;
            register BigReal sortValue1 = list1ptr->sortValue;
            register int index0 = list0ptr->index;
            register int index1 = list1ptr->index;

            // While both lists have at least one element in them, compare the
            //   heads of each list and place the smaller of the two in the
            //   destination array.
            while (list0ptr < list0ptr_end && list1ptr < list1ptr_end) {

	      // Compare the values
              register bool test = (sortValue0 < sortValue1);

              // Place the "winner" in the destination array
              dstptr->sortValue = TERNARY_ASSIGN(test, sortValue0, sortValue1);
              dstptr->index = TERNARY_ASSIGN(test, index0, index1);
              dstptr++;

              // Update the pointers
              list0ptr += TERNARY_ASSIGN(test, 1, 0);
              list1ptr += TERNARY_ASSIGN(test, 0, 1);

              // Refill the sortValue and index register
              // NOTE: These memory locations are likely to be in cache
              sortValue0 = list0ptr->sortValue;
              sortValue1 = list1ptr->sortValue;
              index0 = list0ptr->index;
              index1 = list1ptr->index;

	    } // end while (list0ptr < list0ptr_end && list1ptr < list1ptr_end)

            // NOTE: At this point, at least one of the lists is empty so no
            //   more than one of the loops will be executed.

            // Drain the remaining elements from list0
            while (list0ptr < list0ptr_end) {

	      // Place the value into the destination array
	      dstptr->sortValue = sortValue0;
              dstptr->index = index0;
              dstptr++;

              // Load the next entry in list0
              list0ptr++;
              sortValue0 = list0ptr->sortValue;
              index0 = list0ptr->index;

	    } // end while (list0ptr < list0ptr_end)

            // Drain the remaining elements from list1
            while (list1ptr < list1ptr_end) {

	      // Place the value into the destination array
	      dstptr->sortValue = sortValue1;
              dstptr->index = index1;
              dstptr++;

              // Load the next entry in list1
              list1ptr++;
              sortValue1 = list1ptr->sortValue;
              index1 = list1ptr->index;

	    } // end while (list1ptr < list1ptr_end)

	  } // end while (firstListOffset < p_1_sortValues_len) {

          // Swap the dstArray and srcArray pointers
          register SortEntry* tmpPtr = dstArray;
          dstArray = srcArray;
          srcArray = tmpPtr;

          // Double the subListSize
          subListSize <<= 1;

	}  // end while (subListSize < p_1_sortValues_len)

        // Set the sort values pointers (NOTE: srcArray and dstArray are
        //   swapped at the end of each iteration of the merge sort outer-loop).
        p_1_sortValues_fixg = dstArray;
        p_1_sortValues = srcArray;

        #endif

      #endif
    }


    // Atom list 0
    {
      // Loop through the first list of atoms (p[0])
      // NOTE: This list will NOT be sorted, so the indexes into p_0
      //   and p_0_sortValues match up.  Also, all atoms in this list
      //   (the outer-loop list, need a sort value so do all atoms).

      register int i = 0;
      register unsigned int ngia = p_0->nonbondedGroupIsAtom;
      register unsigned int hgs = p_0->hydrogenGroupSize;
      register BigReal p_x = p_0->position.x;
      register BigReal p_y = p_0->position.y;
      register BigReal p_z = p_0->position.z;

      // Try to take advantage of multiply-add instructions by pulling out the
      //   'V_length2 * V_length_multiplier' term in the sortVal calculation.
      //   So...  '(a2 - b2 + V_length2) * V_length_multiplier' becomes...
      //          '((a2 - b2) * V_length_multiplier) + (V_length2 * V_length_multiplier)'
      register const BigReal sortVal_addAmmount = V_length2 * V_length_multiplier;

      while (i < i_upper) {

	// Advance i... NOTE: ngia is either 0 or 1 so if ngia is set '1' is
	//   added to 'i'.  Otherwise, the value of 'hgs' is added to 'i'.
	register BigReal* p_0_sortValStorePtr = p_0_sortValues + i;
	i += (ngia) + ((1 - ngia) * hgs);

	// Update p_i_next to point to the atom for the next iteration and begin
	//   loading the 'ngia' and 'hgs' values for that atom.
        register const CompAtom* p_i_next = p_0 + i;
        ngia = p_i_next->nonbondedGroupIsAtom;
        hgs = p_i_next->hydrogenGroupSize;

        // Use the position values then start loading the next atom's position values.
        register BigReal a_x_diff = (p_0_avgX - p_x);
        register BigReal a_y_diff = (p_0_avgY - p_y);
        register BigReal a_z_diff = (p_0_avgZ - p_z);
        register BigReal b_x_diff = (p_1_avgX - p_x);
        register BigReal b_y_diff = (p_1_avgY - p_y);
        register BigReal b_z_diff = (p_1_avgZ - p_z);
        p_x = p_i_next->position.x;
        p_y = p_i_next->position.y;
        p_z = p_i_next->position.z;

        // Calculate and store the sort value for this atom (NOTE: c = V_length and c^2 = V_length2)
        register BigReal a2 = (a_x_diff * a_x_diff) + (a_y_diff * a_y_diff) + (a_z_diff * a_z_diff);
        register BigReal b2 = (b_x_diff * b_x_diff) + (b_y_diff * b_y_diff) + (b_z_diff * b_z_diff);
        register BigReal sortVal = (a2 - b2) * V_length_multiplier + sortVal_addAmmount;
        *p_0_sortValStorePtr = sortVal;
      }

      p_0_sortValues_len = i_upper;
    }

  }  // end if (savePairlists || !usePairlists)
  #endif


  // DMK - Atom Sort
  // NOTE: These arrays aren't used for pair computes that are spacially
  //   sorting atoms.
  #if ! (NAMD_ComputeNonbonded_SortAtoms != 0 && ( 0 PAIR( + 1 ) ) )
    NBWORKARRAY(plint,grouplist,arraysize);
    NBWORKARRAY(plint,fixglist,arraysize);
  #endif

  NBWORKARRAY(plint,goodglist,arraysize);
  NBWORKARRAY(plint,pairlistx,arraysize);
  NBWORKARRAY(plint,pairlistm,arraysize);
  NBWORKARRAY(plint,pairlist,arraysize);
  NBWORKARRAY(plint,pairlist2,arraysize);
  ALCH(
  NBWORKARRAY(plint,pairlistnAlch,arraysize);
  NBWORKARRAY(plint,pairlistnA0,arraysize);
  NBWORKARRAY(plint,pairlistnA1,arraysize);
  NBWORKARRAY(plint,pairlistxA1,arraysize);
  NBWORKARRAY(plint,pairlistmA1,arraysize);
  NBWORKARRAY(plint,pairlistnA2,arraysize);
  NBWORKARRAY(plint,pairlistxA2,arraysize);
  NBWORKARRAY(plint,pairlistmA2,arraysize);
  )

  NBWORKARRAY(short,vdwtype_array,j_upper+5);

  
#ifdef MEM_OPT_VERSION
  for (j = 0; j < j_upper; ++j){
    vdwtype_array[j] = pExt_1[j].vdwType;
  }
#else
  const Atom *atomlist = mol->getAtoms();
#ifdef ARCH_POWERPC
#pragma disjoint (*atomlist, *vdwtype_array)
#pragma disjoint (*p_1, *vdwtype_array)
#pragma unroll(4)
#endif
  for (j = 0; j < j_upper; ++j) {
    int id = p_1[j].id;
    vdwtype_array [j] = atomlist[id].vdw_type;
  }
#endif

  int fixg_upper = 0;
  int g_upper = 0;

  if ( savePairlists || ! usePairlists ) {


    // DMK - Atom Sort
    #if NAMD_ComputeNonbonded_SortAtoms != 0 && ( 0 PAIR ( + 1 ) )

      // Create a sorted list of non-fixed groups
      register int fixg = 0;
      for (int tmpI = 0; tmpI < p_1_sortValues_len; tmpI++) {
        register SortEntry* p_1_sortEntry = p_1_sortValues + tmpI;
        register int p_1_index = p_1_sortEntry->index;
        if (!p_1[p_1_index].groupFixed) {
          register SortEntry* p_1_sortEntry_fixg = p_1_sortValues_fixg + p_1_sortValues_fixg_len;
          p_1_sortEntry_fixg->index = p_1_sortEntry->index;
          p_1_sortEntry_fixg->sortValue = p_1_sortEntry->sortValue;
          p_1_sortValues_fixg_len++;
	}
      }

    #else

      register int g = 0;
      for ( j = 0; j < j_upper; ++j ) {
        if ( p_1[j].hydrogenGroupSize || p_1[j].nonbondedGroupIsAtom ) {
          grouplist[g++] = j;
        }
      }
      g_upper = g;
      if ( g_upper ) grouplist[g_upper] = grouplist[g_upper-1];
      int fixg = 0;

      if ( fixedAtomsOn ) {
        for ( g = 0; g < g_upper; ++g ) {
          j = grouplist[g];
          if ( ! p_1[j].groupFixed ) {
            fixglist[fixg++] = j;
          }
        }
      }

      fixg_upper = fixg;
      if ( fixg_upper ) fixglist[fixg_upper] = fixglist[fixg_upper-1];

    #endif // NAMD_ComputeNonbonded_AtomSort != 0



    *(pairlists.newlist(1)) = i_upper;
    pairlists.newsize(1);

  } else { // if ( savePairlists || ! usePairlists )

    plint *i_upper_check;
    int i_upper_check_count;
    pairlists.nextlist(&i_upper_check,&i_upper_check_count);
    if ( i_upper_check[0] != i_upper )
      NAMD_bug("pairlist i_upper mismatch!");

  } // if ( savePairlists || ! usePairlists )

  SELF(
  int j_hgroup = 0;
  int g_lower = 0;
  int fixg_lower = 0;
  )
  int pairlistindex=0;
  int pairlistoffset=0;

  

#if ( SHORT( FAST( 1+ ) ) 0 )
#if ( PAIR( 1+ ) 0 )
    Force *f_0 = params->ff[0];
    Force *f_1 = params->ff[1];
#else
#define f_1 f_0
    NBWORKARRAY(Force,f_0,i_upper)
    memset( (void*) f_0, 0, i_upper * sizeof(Force) );
#endif
#endif
#if ( FULL( 1+ ) 0 )
#if ( PAIR( 1+ ) 0 )
    Force *fullf_0 = params->fullf[0];
    Force *fullf_1 = params->fullf[1];
#else
#define fullf_1 fullf_0
    NBWORKARRAY(Force,fullf_0,i_upper);
    memset( (void*) fullf_0, 0, i_upper * sizeof(Force) );
#endif
#endif
    

  int numParts = params->numParts;
  int myPart = params->minPart;
  int groupCount = 0;

  for ( i = 0; i < (i_upper SELF(- 1)); ++i )
  {
    const CompAtom &p_i = p_0[i];
#ifdef MEM_OPT_VERSION
    const CompAtomExt &pExt_i = pExt_0[i];
#endif
    if ( p_i.hydrogenGroupSize ) {
      //save current group count
      int curgrpcount = groupCount;      
      //increment group count
      groupCount ++;

      if (groupCount >= numParts)
	groupCount = 0;
      
      if ( curgrpcount != myPart ) {
        i += p_i.hydrogenGroupSize - 1;
	
	//Power PC alignment constraint
#ifdef ARCH_POWERPC
	__dcbt((void *) &(p_0[i+1]));
#endif
        continue;
      }
    }

    register const BigReal p_i_x = p_i.position.x + offset_x;
    register const BigReal p_i_y = p_i.position.y + offset_y;
    register const BigReal p_i_z = p_i.position.z + offset_z;

    ALCH(const int p_i_partition = p_i.partition;)

    PPROF(
        const int p_i_partition = p_i.partition;
        int n1 = (int)floor((p_i.position.z-pressureProfileMin)*invThickness);
        pp_clamp(n1, pressureProfileSlabs);
        )

  SELF ( if ( p_i.hydrogenGroupSize ) j_hgroup = i + p_i.hydrogenGroupSize; )

  if ( savePairlists || ! usePairlists ) {

    if ( ! savePairlists ) pairlists.reset();  // limit space usage

    #ifdef MEM_OPT_VERSION
    const ExclusionCheck *exclcheck = mol->get_excl_check_for_idx(pExt_i.exclId);        
    const int excl_min = p_i.id + exclcheck->min;
    const int excl_max = p_i.id + exclcheck->max;
    #else
    const ExclusionCheck *exclcheck = mol->get_excl_check_for_atom(p_i.id);
    const int excl_min = exclcheck->min;
    const int excl_max = exclcheck->max;
    #endif
    const char * excl_flags_var;
    if ( exclcheck->flags ) excl_flags_var = exclcheck->flags - excl_min;
    else {  // need to build list on the fly

    //TODO: Should change later!!!!!!!!!! --Chao Mei
    //Now just for the sake of passing compilation
    #ifndef MEM_OPT_VERSION 
      if ( excl_flags_buff ) {
        int nl,l;
        nl = full_excl[0] + 1;
        for ( l=1; l<nl; ++l ) excl_flags_buff[full_excl[l]] = 0;
        nl = mod_excl[0] + 1;
        for ( l=1; l<nl; ++l ) excl_flags_buff[mod_excl[l]] = 0;
      } else {
        excl_flags_buff = new char[mol->numAtoms];
        memset( (void*) excl_flags_buff, 0, mol->numAtoms);
      }
      int nl,l;
      full_excl = mol->get_full_exclusions_for_atom(p_i.id);
      nl = full_excl[0] + 1;
      for ( l=1; l<nl; ++l ) excl_flags_buff[full_excl[l]] = EXCHCK_FULL;
      mod_excl = mol->get_mod_exclusions_for_atom(p_i.id);
      nl = mod_excl[0] + 1;
      for ( l=1; l<nl; ++l ) excl_flags_buff[mod_excl[l]] = EXCHCK_MOD;
      excl_flags_var = excl_flags_buff;
    #endif

    }
    const char * const excl_flags = excl_flags_var;

  if (p_i.hydrogenGroupSize || p_i.nonbondedGroupIsAtom) {

    pairlistindex = 0;	// initialize with 0 elements
    pairlistoffset=0;
    const int groupfixed = ( fixedAtomsOn && p_i.groupFixed );

    // If patch divisions are not made by hydrogen groups, then
    // hydrogenGroupSize is set to 1 for all atoms.  Thus we can
    // carry on as if we did have groups - only less efficiently.
    // An optimization in this case is to not rebuild the temporary
    // pairlist but to include every atom in it.  This should be a
    // a very minor expense.

    SELF
    (
      if ( p_i.hydrogenGroupSize ) {
        // exclude child hydrogens of i
        // j_hgroup = i + p_i.hydrogenGroupSize;  (moved above)
        while ( g_lower < g_upper &&
                grouplist[g_lower] < j_hgroup ) ++g_lower;
        while ( fixg_lower < fixg_upper &&
                fixglist[fixg_lower] < j_hgroup ) ++fixg_lower;
      }
      // add all child or sister hydrogens of i
      for ( j = i + 1; j < j_hgroup; ++j ) {
	pairlist[pairlistindex++] = j;
      }
    )

    // add remaining atoms to pairlist via hydrogen groups
    register plint *pli = pairlist + pairlistindex;

    {
      // DMK - Atom Sort
      // Modify the values of g and gu based on the added information
      //   of the linear projections (sort values) information.
      #if NAMD_ComputeNonbonded_SortAtoms != 0 && ( 0 PAIR ( + 1 ) )

        register SortEntry* sortValues = ( groupfixed ? p_1_sortValues_fixg : p_1_sortValues );
        register int g = 0;
        const int gu = ( groupfixed ? p_1_sortValues_fixg_len : p_1_sortValues_len );

        register BigReal p_i_sortValue = p_0_sortValues[i];
        const BigReal p_i_sortValue_plus_windowRadius = p_i_sortValue + atomSort_windowRadius;

      #else

        register plint *gli = goodglist;
        const plint *glist = ( groupfixed ? fixglist : grouplist );
        SELF( const int gl = ( groupfixed ? fixg_lower : g_lower ); )
        const int gu = ( groupfixed ? fixg_upper : g_upper );
        register int g = PAIR(0) SELF(gl);

      #endif


      if ( g < gu ) {
	int hu = 0;
#if defined(NAMD_SSE) && defined(__INTEL_COMPILER) && defined(__SSE2__)
	if ( gu - g  >  6 ) { 

	  // DMK - Atom Sort
          #if NAMD_ComputeNonbonded_SortAtoms != 0 && ( 0 PAIR ( + 1 ) )

	    register SortEntry* sortEntry0 = sortValues + g;
	    register SortEntry* sortEntry1 = sortValues + g + 1;
            register int jprev0 = sortEntry0->index;
	    register int jprev1 = sortEntry1->index;

            // Test if we need to loop again
            // Only if both of the atoms in the next iteration are
            //   outside of the bounds.  If one is within, it must
            //   be added to goodglist next iteration so keep going.
            bool test2 = ((p_i_sortValue_plus_windowRadius < sortEntry0->sortValue)
                       && (p_i_sortValue_plus_windowRadius < sortEntry1->sortValue));
            g += (test2 * gu);  // Either add '0' or 'gu'

          #else

            register int jprev0 = glist[g    ];
	    register int jprev1 = glist[g + 1];

          #endif
	  
	  register int j0;
	  register int j1;

          __m128d PJ_X_01 = _mm_set_pd(p_1[jprev1].position.x, p_1[jprev0].position.x);
          __m128d PJ_Y_01 = _mm_set_pd(p_1[jprev1].position.y, p_1[jprev0].position.y);
          __m128d PJ_Z_01 = _mm_set_pd(p_1[jprev1].position.z, p_1[jprev0].position.z);

          // these don't change here, so we could move them into outer scope
          const __m128d R2_DELTA = _mm_set1_pd(r2_delta);	 
          const __m128d P_I_X = _mm_set1_pd(p_i_x);
          const __m128d P_I_Y = _mm_set1_pd(p_i_y);
          const __m128d P_I_Z = _mm_set1_pd(p_i_z);
 
	  g += 2;
	  for ( ; g < gu - 2; g +=2 ) {
	    // compute 1d distance, 2-way parallel	 
	    j0     =  jprev0;
	    j1     =  jprev1;

            __m128d T_01 = _mm_sub_pd(P_I_X, PJ_X_01);
            __m128d R2_01 = _mm_add_pd(_mm_mul_pd(T_01, T_01), R2_DELTA);
            T_01 = _mm_sub_pd(P_I_Y, PJ_Y_01);
            R2_01 = _mm_add_pd(R2_01, _mm_mul_pd(T_01, T_01));
            T_01 = _mm_sub_pd(P_I_Z, PJ_Z_01);
            R2_01 = _mm_add_pd(R2_01, _mm_mul_pd(T_01, T_01));
	    
            // DMK - Atom Sort
            #if NAMD_ComputeNonbonded_SortAtoms != 0 && ( 0 PAIR ( + 1 ) )

	      sortEntry0 = sortValues + g;
	      sortEntry1 = sortValues + g + 1;
              jprev0 = sortEntry0->index;
	      jprev1 = sortEntry1->index;

              // Test if we need to loop again
              // Only if both of the atoms in the next iteration are
              //   outside of the bounds.  If one is within, it must
              //   be added to goodglist next iteration so keep going.
              bool test2 = ((p_i_sortValue_plus_windowRadius < sortEntry0->sortValue)
                         && (p_i_sortValue_plus_windowRadius < sortEntry1->sortValue));
              g += (test2 * gu);  // Either add '0' or 'gu'

            #else

	      jprev0     =  glist[g  ];
	      jprev1     =  glist[g+1];

            #endif
	   
            PJ_X_01 = _mm_set_pd(p_1[jprev1].position.x, p_1[jprev0].position.x);
            PJ_Y_01 = _mm_set_pd(p_1[jprev1].position.y, p_1[jprev0].position.y);
            PJ_Z_01 = _mm_set_pd(p_1[jprev1].position.z, p_1[jprev0].position.z);

            __declspec(align(16)) double r2_01[2];
            _mm_store_pd(r2_01, R2_01); // 16-byte-aligned store

            // XXX these could possibly benefit from SSE-based conditionals
	    bool test0 = ( r2_01[0] < groupplcutoff2 );
	    bool test1 = ( r2_01[1] < groupplcutoff2 ); 
	    
	    //removing ifs benefits on many architectures
	    //as the extra stores will only warm the cache up
	    goodglist [ hu         ] = j0;
	    goodglist [ hu + test0 ] = j1;
	    
	    hu += test0 + test1;
	  }
	  g-=2;
	}
#else
	if ( gu - g  >  6 ) { 

	  // DMK - Atom Sort
          #if NAMD_ComputeNonbonded_SortAtoms != 0 && ( 0 PAIR ( + 1 ) )

	    register SortEntry* sortEntry0 = sortValues + g;
	    register SortEntry* sortEntry1 = sortValues + g + 1;
            register int jprev0 = sortEntry0->index;
	    register int jprev1 = sortEntry1->index;

            // Test if we need to loop again
            // Only if both of the atoms in the next iteration are
            //   outside of the bounds.  If one is within, it must
            //   be added to goodglist next iteration so keep going.
            bool test2 = ((p_i_sortValue_plus_windowRadius < sortEntry0->sortValue)
                       && (p_i_sortValue_plus_windowRadius < sortEntry1->sortValue));
            g += (test2 * gu);  // Either add '0' or 'gu'

          #else

            register int jprev0 = glist[g    ];
	    register int jprev1 = glist[g + 1];

          #endif
	  
	  register  int j0; 
	  register  int j1; 
	  
	  register  BigReal pj_x_0, pj_x_1; 
	  register  BigReal pj_y_0, pj_y_1; 
	  register  BigReal pj_z_0, pj_z_1; 
	  register  BigReal t_0, t_1, r2_0, r2_1;
	  
	  pj_x_0 = p_1[jprev0].position.x;
	  pj_x_1 = p_1[jprev1].position.x;  
	  
	  pj_y_0 = p_1[jprev0].position.y; 
	  pj_y_1 = p_1[jprev1].position.y;  
	  
	  pj_z_0 = p_1[jprev0].position.z; 
	  pj_z_1 = p_1[jprev1].position.z;
	  
	  g += 2;
	  for ( ; g < gu - 2; g +=2 ) {
	    // compute 1d distance, 2-way parallel	 
	    j0     =  jprev0;
	    j1     =  jprev1;
	    
	    t_0    =  p_i_x - pj_x_0;
	    t_1    =  p_i_x - pj_x_1;
	    r2_0   =  t_0 * t_0;
	    r2_1   =  t_1 * t_1;
	    
	    t_0    =  p_i_y - pj_y_0;
	    t_1    =  p_i_y - pj_y_1;
	    r2_0  +=  t_0 * t_0;
	    r2_1  +=  t_1 * t_1;
	    
	    t_0    =  p_i_z - pj_z_0;
	    t_1    =  p_i_z - pj_z_1;
	    r2_0  +=  t_0 * t_0;
	    r2_1  +=  t_1 * t_1;
	    
            // DMK - Atom Sort
            #if NAMD_ComputeNonbonded_SortAtoms != 0 && ( 0 PAIR ( + 1 ) )

	      sortEntry0 = sortValues + g;
	      sortEntry1 = sortValues + g + 1;
              jprev0 = sortEntry0->index;
	      jprev1 = sortEntry1->index;

              // Test if we need to loop again
              // Only if both of the atoms in the next iteration are
              //   outside of the bounds.  If one is within, it must
              //   be added to goodglist next iteration so keep going.
              bool test2 = ((p_i_sortValue_plus_windowRadius < sortEntry0->sortValue)
                         && (p_i_sortValue_plus_windowRadius < sortEntry1->sortValue));
              g += (test2 * gu);  // Either add '0' or 'gu'

            #else

	      jprev0     =  glist[g  ];
	      jprev1     =  glist[g+1];

            #endif
	    
	    pj_x_0     =  p_1[jprev0].position.x;
	    pj_x_1     =  p_1[jprev1].position.x;
	    pj_y_0     =  p_1[jprev0].position.y; 
	    pj_y_1     =  p_1[jprev1].position.y;
	    pj_z_0     =  p_1[jprev0].position.z; 
	    pj_z_1     =  p_1[jprev1].position.z;
	    
	    bool test0 = ( r2_0 < groupplcutoff2 );
	    bool test1 = ( r2_1 < groupplcutoff2 ); 
	    
	    //removing ifs benefits on many architectures
	    //as the extra stores will only warm the cache up
	    goodglist [ hu         ] = j0;
	    goodglist [ hu + test0 ] = j1;
	    
	    hu += test0 + test1;
	  }
	  g-=2;
	}
#endif
	
	for (; g < gu; g++) {

          // DMK - Atom Sort
          #if NAMD_ComputeNonbonded_SortAtoms != 0 && ( 0 PAIR ( + 1 ) )

	    register SortEntry* sortEntry = sortValues + g;
            register int j = sortEntry->index;

            // Test if we need to loop again
            // NOTE: To avoid a branch, always do this iteration's distance
            //   calculation and just test this sort value to see if the
            //   rest of the iterations can be skiped (vs. a 'break' if
            //   this one can be skiped).
            bool test = (p_i_sortValue_plus_windowRadius < sortEntry->sortValue);
            g += (test * gu);  // Either add '0' or 'gu'

          #else

	    int j = glist[g];

          #endif

	  BigReal p_j_x = p_1[j].position.x;
	  BigReal p_j_y = p_1[j].position.y;
	  BigReal p_j_z = p_1[j].position.z;
	  
	  BigReal r2 = p_i_x - p_j_x;
	  r2 *= r2;
	  BigReal t2 = p_i_y - p_j_y;
	  r2 += t2 * t2;
	  t2 = p_i_z - p_j_z;
	  r2 += t2 * t2;

	  if ( r2 <= groupplcutoff2 ) 
	    goodglist[hu ++] = j; 
	}

	for ( int h=0; h<hu; ++h ) {
          int j = goodglist[h];
          int hgs = ( p_1[j].nonbondedGroupIsAtom ? 1 :
		      p_1[j].hydrogenGroupSize );
	  pli[0] = j;   // copy over the next four in any case
	  pli[1] = j+1;
	  pli[2] = j+2;
	  pli[3] = j+3; // assume hgs <= 4
          pli += hgs;
	}
      }
    }

    pairlistindex = pli - pairlist;
    // make sure padded element on pairlist points to real data
    if ( pairlistindex ) {
       pairlist[pairlistindex] = pairlist[pairlistindex-1];
    } /* PAIR( else {  // skip empty loops if no pairs were found
       int hgs = ( p_i.nonbondedGroupIsAtom ? 1 : p_i.hydrogenGroupSize );
       i += hgs - 1;
       continue;
    } ) */
  } // if i is hydrogen group parent
  SELF
    (
    // self-comparisions require list to be incremented
    // pair-comparisions use entire list (pairlistoffset is 0)
    else pairlistoffset++;
    )

    const int atomfixed = ( fixedAtomsOn && p_i.atomFixed );

    register plint *pli = pairlist2;
#if 1 ALCH(-1)
    plint *pairlistn = pairlists.newlist(j_upper + 5 + 1 + 5) SELF(+ 1);
#else
    plint *pairlistn = pairlistnAlch;
#endif
    SELF( plint &pairlistn_skip = *(pairlistn-1); )
    register plint *plin = pairlistn;


    INT(
    if ( pairInteractionOn ) {
      const int ifep_type = p_i.partition;
      if (pairInteractionSelf) {
        if (ifep_type != 1) continue;
        for (int k=pairlistoffset; k<pairlistindex; k++) {
          j = pairlist[k];
          const int jfep_type = p_1[j].partition;
          // for pair-self, both atoms must be in group 1.
          if (jfep_type == 1) {
            *(pli++) = j;
          }
        }
      } else {
        if (ifep_type != 1 && ifep_type != 2) continue;
        for (int k=pairlistoffset; k<pairlistindex; k++) {
          j = pairlist[k];
          const int jfep_type = p_1[j].partition;
          // for pair, must have one from each group.
          if (ifep_type + jfep_type == 3) {
            *(pli++) = j;
          }
        }
      }
      int npair2_int = pli - pairlist2;
      pli = pairlist2;
      for (int k=0; k<npair2_int; k++) {
        j = pairlist2[k];
        BigReal p_j_x = p_1[j].position.x;
	BigReal r2 = p_i_x - p_j_x;
	r2 *= r2;
        BigReal p_j_y = p_1[j].position.y;
	BigReal t2 = p_i_y - p_j_y;
	r2 += t2 * t2;
        BigReal p_j_z = p_1[j].position.z;
	t2 = p_i_z - p_j_z;
	r2 += t2 * t2;
	if ( ( ! (atomfixed && p_1[j].atomFixed) ) && (r2 <= plcutoff2) ) {
          int atom2 = p_1[j].id;
          if ( atom2 >= excl_min && atom2 <= excl_max ) *(pli++) = j;
          else *(plin++) = j;
        }
      }
    } else
    )
    if ( atomfixed ) {
      for (int k=pairlistoffset; k<pairlistindex; k++) {
        j = pairlist[k];
        BigReal p_j_x = p_1[j].position.x;
	BigReal r2 = p_i_x - p_j_x;
	r2 *= r2;
        BigReal p_j_y = p_1[j].position.y;
	BigReal t2 = p_i_y - p_j_y;
	r2 += t2 * t2;
        BigReal p_j_z = p_1[j].position.z;
	t2 = p_i_z - p_j_z;
	r2 += t2 * t2;
	if ( (! p_1[j].atomFixed) && (r2 <= plcutoff2) ) {
          int atom2 = p_1[j].id;
          if ( atom2 >= excl_min && atom2 <= excl_max ) *(pli++) = j;
          else *(plin++) = j;
        }
      }
    } else {
      int k = pairlistoffset;
      int ku = pairlistindex;
      if ( k < ku ) {
#if defined(NAMD_SSE) && defined(__INTEL_COMPILER) && defined(__SSE2__)
	if ( ku - k  >  6 ) { 	   
	  register  int jprev0 = pairlist [k    ];
	  register  int jprev1 = pairlist [k + 1];
	  
	  register  int j0; 
	  register  int j1; 

          __m128d PJ_X_01 = _mm_set_pd(p_1[jprev1].position.x, p_1[jprev0].position.x);
          __m128d PJ_Y_01 = _mm_set_pd(p_1[jprev1].position.y, p_1[jprev0].position.y);
          __m128d PJ_Z_01 = _mm_set_pd(p_1[jprev1].position.z, p_1[jprev0].position.z);

          // these don't change here, so we could move them into outer scope
          const __m128d R2_DELTA = _mm_set1_pd(r2_delta);
          const __m128d P_I_X = _mm_set1_pd(p_i_x);
          const __m128d P_I_Y = _mm_set1_pd(p_i_y);
          const __m128d P_I_Z = _mm_set1_pd(p_i_z);
	  
	  int atom2_0 = p_1[jprev0].id;
	  int atom2_1 = p_1[jprev1].id;
	  
	  k += 2;
	  for ( ; k < ku - 2; k +=2 ) {
	    // compute 1d distance, 2-way parallel	 
	    j0     =  jprev0;
	    j1     =  jprev1;
	    
            __m128d T_01 = _mm_sub_pd(P_I_X, PJ_X_01);
            __m128d R2_01 = _mm_add_pd(_mm_mul_pd(T_01, T_01), R2_DELTA);
            T_01 = _mm_sub_pd(P_I_Y, PJ_Y_01);
            R2_01 = _mm_add_pd(R2_01, _mm_mul_pd(T_01, T_01));
            T_01 = _mm_sub_pd(P_I_Z, PJ_Z_01);
            R2_01 = _mm_add_pd(R2_01, _mm_mul_pd(T_01, T_01));
	    
	    jprev0     =  pairlist[k];
	    jprev1     =  pairlist[k+1];
	    
            PJ_X_01 = _mm_set_pd(p_1[jprev1].position.x, p_1[jprev0].position.x);
            PJ_Y_01 = _mm_set_pd(p_1[jprev1].position.y, p_1[jprev0].position.y);
            PJ_Z_01 = _mm_set_pd(p_1[jprev1].position.z, p_1[jprev0].position.z);

            __declspec(align(16)) double r2_01[2];
            _mm_store_pd(r2_01, R2_01); // 16-byte-aligned store
	    
	    if (r2_01[0] <= plcutoff2) {
	      if ( atom2_0 >= excl_min && atom2_0 <= excl_max ) 
		*(pli++) = j0;
	      else 
		*(plin++) = j0;
	    }
	    atom2_0 = p_1[jprev0].id;
	    
	    if (r2_01[1] <= plcutoff2) {
	      if ( atom2_1 >= excl_min && atom2_1 <= excl_max ) 
		*(pli++) = j1;
	      else 
		*(plin++) = j1;
	    }
	    atom2_1 = p_1[jprev1].id;	    
	  }
	  k-=2;
	}       
#else
	if ( ku - k  >  6 ) { 	   
	  register  int jprev0 = pairlist [k];
	  register  int jprev1 = pairlist [k + 1];
	  
	  register  int j0; 
	  register  int j1; 
	  
	  register  BigReal pj_x_0, pj_x_1; 
	  register  BigReal pj_y_0, pj_y_1; 
	  register  BigReal pj_z_0, pj_z_1; 
	  register  BigReal t_0, t_1, r2_0, r2_1;
	  
	  pj_x_0 = p_1[jprev0].position.x;
	  pj_x_1 = p_1[jprev1].position.x;  
	  
	  pj_y_0 = p_1[jprev0].position.y; 
	  pj_y_1 = p_1[jprev1].position.y;  
	  
	  pj_z_0 = p_1[jprev0].position.z; 
	  pj_z_1 = p_1[jprev1].position.z;
	  
	  int atom2_0 = p_1[jprev0].id;
	  int atom2_1 = p_1[jprev1].id;
	  
	  k += 2;
	  for ( ; k < ku - 2; k +=2 ) {
	    // compute 1d distance, 2-way parallel	 
	    j0     =  jprev0;
	    j1     =  jprev1;
	    
	    t_0    =  p_i_x - pj_x_0;
	    t_1    =  p_i_x - pj_x_1;
	    r2_0   =  t_0 * t_0;
	    r2_1   =  t_1 * t_1;
	    
	    t_0    =  p_i_y - pj_y_0;
	    t_1    =  p_i_y - pj_y_1;
	    r2_0  +=  t_0 * t_0;
	    r2_1  +=  t_1 * t_1;
	    
	    t_0    =  p_i_z - pj_z_0;
	    t_1    =  p_i_z - pj_z_1;
	    r2_0  +=  t_0 * t_0;
	    r2_1  +=  t_1 * t_1;
	    
	    jprev0     =  pairlist[k];
	    jprev1     =  pairlist[k+1];
	    
	    pj_x_0     =  p_1[jprev0].position.x;
	    pj_x_1     =  p_1[jprev1].position.x;
	    pj_y_0     =  p_1[jprev0].position.y; 
	    pj_y_1     =  p_1[jprev1].position.y;
	    pj_z_0     =  p_1[jprev0].position.z; 
	    pj_z_1     =  p_1[jprev1].position.z;
	    
	    if (r2_0 <= plcutoff2) {
	      if ( atom2_0 >= excl_min && atom2_0 <= excl_max ) 
		*(pli++) = j0;
	      else 
		*(plin++) = j0;
	    }
	    atom2_0 = p_1[jprev0].id;
	    
	    if (r2_1 <= plcutoff2) {
	      if ( atom2_1 >= excl_min && atom2_1 <= excl_max ) 
		*(pli++) = j1;
	      else 
		*(plin++) = j1;
	     }
	    atom2_1 = p_1[jprev1].id;	    
	  }
	  k-=2;
	}       
#endif

	for (; k < ku; k++) {
	  int j = pairlist[k];
	  int atom2 = p_1[j].id;
	  
	  BigReal p_j_x = p_1[j].position.x;
	  BigReal p_j_y = p_1[j].position.y;
	  BigReal p_j_z = p_1[j].position.z;
	  
	  BigReal r2 = p_i_x - p_j_x;
	  r2 *= r2;
	  BigReal t2 = p_i_y - p_j_y;
	  r2 += t2 * t2;
	  t2 = p_i_z - p_j_z;
	  r2 += t2 * t2;
	  
	  if (r2 <= plcutoff2) {
	    if ( atom2 >= excl_min && atom2 <= excl_max ) 
	      *(pli++) = j;
	    else 
	      *(plin++) = j;
	  }
	}
      }
    }
    int npair2 = pli - pairlist2;
    if ( npair2 ) pairlist2[npair2] = pairlist2[npair2-1];

    plint *plix = pairlistx;
    plint *plim = pairlistm;
    plint *pln = pairlistn;
    ALCH(
    plint *plinA1 = pairlistnA1;
    plint *plixA1 = pairlistxA1;
    plint *plimA1 = pairlistmA1;
    plint *plinA2 = pairlistnA2;
    plint *plixA2 = pairlistxA2;
    plint *plimA2 = pairlistmA2;
    )
    int k=0;
#if 1 ALCH(-1)
    SELF(
    for (; pln < plin && *pln < j_hgroup; ++pln) {
      *(plix++) = *pln;  // --exclChecksum;
    }
    pairlistn_skip = pln - pairlistn;
    for (; k < npair2 && pairlist2[k] < j_hgroup; ++k) {
      *(plix++) = pairlist2[k];  // --exclChecksum;
    }
    )
#endif
ALCH(
    SELF(
    for (; pln < plin && *pln < j_hgroup; ++pln) {
      switch (pswitchTable[4*p_i_partition]) { //p_i_partition + 3*p_1[i].partition
      case 0: *(plix++) = *pln;  break;
      case 1: *(plixA1++) = *pln; break;
      case 2: *(plixA2++) = *pln; break;
      default: NAMD_die("Alchemical pairlist error"); break;
      }
    }
    for (; k < npair2 && pairlist2[k] < j_hgroup; ++k) {
      switch (pswitchTable[4*p_i_partition]) {
      case 0: *(plix++) = pairlist2[k];  break;
      case 1: *(plixA1++) = pairlist2[k]; break;
      case 2: *(plixA2++) = pairlist2[k]; break;
      default: NAMD_die("Alchemical pairlist error"); break;
      }
    }
    )
)
    for (; k < npair2; ++k ) {
      int j = pairlist2[k];
      int atom2 = p_1[j].id;
      int excl_flag = excl_flags[atom2];
      ALCH(int pswitch = pswitchTable[p_i_partition + 3*(p_1[j].partition)];)
      switch ( excl_flag ALCH( + 3 * pswitch)) {
      case 0:  *(plin++) = j;  break;
      case 1:  *(plix++) = j;  break;
      case 2:  *(plim++) = j;  break;
      ALCH(
      case 3:  *(plinA1++) = j; break;
      case 6:  *(plinA2++) = j; break;
      case 5:  *(plimA1++) = j; break;
      case 8:  *(plimA2++) = j; break;
      case 4:  *(plixA1++) = j; break;
      case 7:  *(plixA2++) = j; break;
      )
      }
    }
    // exclChecksum += (plix - pairlistx);
    // exclChecksum += (plim - pairlistm);


#if 1 ALCH(-1)
    npairn = plin - pln;
    pairlistn_save = pln;
    pairlistn_save[npairn] = npairn ? pairlistn_save[npairn-1] : -1;
    pairlists.newsize(plin - pairlistn SELF(+ 1) + 1);
#else
    plint *plinA0 = pairlistnA0;
    int unsortedNpairn = plin - pln;
    for ( k=0; k<unsortedNpairn; ++k ) {
      int j = pln[k];
      switch(pswitchTable[p_i_partition + 3*(p_1[j].partition)]) {
        case 0:  *(plinA0++) = j; break;
        case 1:  *(plinA1++) = j; break;
        case 2:  *(plinA2++) = j; break;
      }
    }
    
    npairn = plinA0 - pairlistnA0;
    // FB preallocation (incl extra for overhead) seems to be necessary
    pairlistn_save = pairlists.newlist(j_upper + 30);
    for ( k=0; k<npairn; ++k ) {
      pairlistn_save[k] = pairlistnA0[k];
    }
    pairlistn_save[k] = k ? pairlistn_save[k-1] : -1;
    pairlists.newsize(npairn + 1);

#endif
    
    npairx = plix - pairlistx;
    pairlistx_save = pairlists.newlist(npairx + 1);
    for ( k=0; k<npairx; ++k ) {
      pairlistx_save[k] = pairlistx[k];
    }
    pairlistx_save[k] = k ? pairlistx_save[k-1] : -1;
    pairlists.newsize(npairx + 1);

    npairm = plim - pairlistm;
    pairlistm_save = pairlists.newlist(npairm + 1);
    for ( k=0; k<npairm; ++k ) {
      pairlistm_save[k] = pairlistm[k];
    }
    pairlistm_save[k] = k ? pairlistm_save[k-1] : -1;
    pairlists.newsize(npairm + 1);

    
#if 0 ALCH(+1)
#define PAIRLISTFROMARRAY(NPAIRS,PL1,PL2,PLSAVE) \
{ \
  (NPAIRS) = (PL2) - (PL1); \
  (PLSAVE) = pairlists.newlist((NPAIRS) + 1); \
  for ( k=0; k<(NPAIRS); ++k ) { \
    (PLSAVE)[k] = (PL1)[k]; \
  } \
  (PLSAVE)[k] = k ? (PLSAVE)[k-1] : -1; \
  pairlists.newsize((NPAIRS) + 1); \
}
  
    PAIRLISTFROMARRAY(npairnA1,pairlistnA1,plinA1,pairlistnA1_save);
    PAIRLISTFROMARRAY(npairxA1,pairlistxA1,plixA1,pairlistxA1_save);
    PAIRLISTFROMARRAY(npairmA1,pairlistmA1,plimA1,pairlistmA1_save);
    PAIRLISTFROMARRAY(npairnA2,pairlistnA2,plinA2,pairlistnA2_save);
    PAIRLISTFROMARRAY(npairxA2,pairlistxA2,plixA2,pairlistxA2_save);
    PAIRLISTFROMARRAY(npairmA2,pairlistmA2,plimA2,pairlistmA2_save);
#undef PAIRLISTFROMARRAY

#endif
    
    
	// PAIR( iout << i << " " << i_upper << " save\n" << endi;)
  } else { // if ( savePairlists || ! usePairlists )
	// PAIR( iout << i << " " << i_upper << " use\n" << endi;)

    pairlists.nextlist(&pairlistn_save,&npairn);  --npairn;
    //if ( npairn > 1000 )
//	iout << i << " " << i_upper << " " << npairn << " n\n" << endi;
#if 1 ALCH(-1)
    SELF(
    int pairlistn_skip = *pairlistn_save;
    pairlistn_save += (pairlistn_skip + 1);
    npairn -= (pairlistn_skip + 1);
    )
#endif

    pairlists.nextlist(&pairlistx_save,&npairx);  --npairx;

    //if ( npairx > 1000 )
//	iout << i << " " << i_upper << " " << npairx << " x\n" << endi;
    // exclChecksum += npairx;
    pairlists.nextlist(&pairlistm_save,&npairm);  --npairm;
    ALCH(
    pairlists.nextlist(&pairlistnA1_save,&npairnA1);  --npairnA1;
    pairlists.nextlist(&pairlistxA1_save,&npairxA1);  --npairxA1;
    pairlists.nextlist(&pairlistmA1_save,&npairmA1);  --npairmA1;
    pairlists.nextlist(&pairlistnA2_save,&npairnA2);  --npairnA2;
    pairlists.nextlist(&pairlistxA2_save,&npairxA2);  --npairxA2;
    pairlists.nextlist(&pairlistmA2_save,&npairmA2);  --npairmA2;
    )
    //if ( npairm > 1000 )
//	iout << i << " " << i_upper << " " << npairm << " m\n" << endi;
    // exclChecksum += npairm;

  } // if ( savePairlists || ! usePairlists )

    LES( BigReal *lambda_table_i =
			lambda_table + (lesFactor+1) * p_i.partition; )

    INT(
    const BigReal force_sign = (
      ( pairInteractionOn && ! pairInteractionSelf ) ?
        ( ( p_i.partition == 1 ) ? 1. : -1. ) : 0.
    );
      
    )

    const BigReal kq_i = COLOUMB * p_i.charge * scaling * dielectric_1;
#ifdef MEM_OPT_VERSION
    const LJTable::TableEntry * const lj_row =
		ljTable->table_row(pExt_i.vdwType);
#else
    const LJTable::TableEntry * const lj_row =
		ljTable->table_row(mol->atomvdwtype(p_i.id));
#endif

    SHORT( FAST( BigReal f_i_x = 0.; ) )
    SHORT( FAST( BigReal f_i_y = 0.; ) )
    SHORT( FAST( BigReal f_i_z = 0.; ) )
    FULL( BigReal fullf_i_x = 0.; )
    FULL( BigReal fullf_i_y = 0.; )
    FULL( BigReal fullf_i_z = 0.; )

    int npairi;
    int k;

    npairi = pairlist_from_pairlist(ComputeNonbondedUtil::cutoff2,
	p_i_x, p_i_y, p_i_z, p_1, pairlistn_save, npairn, pairlisti,
	r2_delta, r2list);

#define NORMAL(X) X
#define EXCLUDED(X)
#define MODIFIED(X)
#include  "ComputeNonbondedBase2.h"
#undef NORMAL
#undef EXCLUDED
#undef MODIFIED

    npairi = pairlist_from_pairlist(ComputeNonbondedUtil::cutoff2,
	p_i_x, p_i_y, p_i_z, p_1, pairlistm_save, npairm, pairlisti,
	r2_delta, r2list);
    exclChecksum += npairi;

#define NORMAL(X)
#define EXCLUDED(X)
#define MODIFIED(X) X
#include  "ComputeNonbondedBase2.h"
#undef NORMAL
#undef EXCLUDED
#undef MODIFIED

#ifdef FULLELECT
    npairi = pairlist_from_pairlist(ComputeNonbondedUtil::cutoff2,
	p_i_x, p_i_y, p_i_z, p_1, pairlistx_save, npairx, pairlisti,
	r2_delta, r2list);
    exclChecksum += npairi;
    SELF(
    for (k=0; k<npairi && pairlisti[k] < j_hgroup; ++k) --exclChecksum;
    )

#undef FAST
#define FAST(X)
#define NORMAL(X)
#define EXCLUDED(X) X
#define MODIFIED(X)
#include  "ComputeNonbondedBase2.h"
#undef FAST
#ifdef SLOWONLY
  #define FAST(X)
#else
  #define FAST(X) X
#endif
#undef NORMAL
#undef EXCLUDED
#undef MODIFIED
#else
    exclChecksum += npairx;
    SELF(
      for (k=0; k<npairx && pairlistx_save[k] < j_hgroup; ++k) --exclChecksum;
    )
#endif

#if 0 ALCH(+1)   // nonbondedbase2 for alchemical-separeted pairlists

    #undef ALCHPAIR
    #define ALCHPAIR(X) X
    #undef NOT_ALCHPAIR
    #define NOT_ALCHPAIR(X)
    BigReal myLambda; FEP(BigReal myLambda2;)
    BigReal myElecLambda;  FEP(BigReal myElecLambda2;)
    BigReal myVdwLambda; FEP(BigReal myVdwLambda2;)
    BigReal myVdwShift; FEP(BigReal myVdwShift2;)
    BigReal alch_vdw_energy; BigReal alch_vdw_force;
    FEP(BigReal alch_vdw_energy_2;) TI(BigReal alch_vdw_dUdl;)
    BigReal shiftedElec; BigReal shiftedElecForce;
    
    /********************************************************************/
    /*******NONBONDEDBASE2 FOR NORMAL INTERACTIONS SCALED BY LAMBDA******/
    /********************************************************************/
    #define NORMAL(X) X
    #define EXCLUDED(X)
    #define MODIFIED(X)

    #define ALCH1(X) X
    #define ALCH2(X)
    npairi = pairlist_from_pairlist(ComputeNonbondedUtil::cutoff2,
            p_i_x, p_i_y, p_i_z, p_1, pairlistnA1_save, npairnA1, pairlisti,
            r2_delta, r2list);
    #include  "ComputeNonbondedBase2.h" // normal, direction 'up'
    #undef ALCH1
    #undef ALCH2

    #define ALCH1(X)
    #define ALCH2(X) X
    npairi = pairlist_from_pairlist(ComputeNonbondedUtil::cutoff2,
            p_i_x, p_i_y, p_i_z, p_1, pairlistnA2_save, npairnA2, pairlisti,
            r2_delta, r2list);
    #include  "ComputeNonbondedBase2.h" // normal, direction 'down'
    #undef ALCH1
    #undef ALCH2

    #undef NORMAL
    #undef EXCLUDED
    #undef MODIFIED
    
    /********************************************************************/
    /******NONBONDEDBASE2 FOR MODIFIED INTERACTIONS SCALED BY LAMBDA*****/
    /********************************************************************/
    #define NORMAL(X)
    #define EXCLUDED(X)
    #define MODIFIED(X) X

    #define ALCH1(X) X
    #define ALCH2(X)
    npairi = pairlist_from_pairlist(ComputeNonbondedUtil::cutoff2,
            p_i_x, p_i_y, p_i_z, p_1, pairlistmA1_save, npairmA1, pairlisti,
            r2_delta, r2list);
        exclChecksum += npairi;
    #include  "ComputeNonbondedBase2.h" // modified, direction 'up'
    #undef ALCH1
    #undef ALCH2

    #define ALCH1(X)
    #define ALCH2(X) X
    npairi = pairlist_from_pairlist(ComputeNonbondedUtil::cutoff2,
            p_i_x, p_i_y, p_i_z, p_1, pairlistmA2_save, npairmA2, pairlisti,
            r2_delta, r2list);
        exclChecksum += npairi;
    #include  "ComputeNonbondedBase2.h" // modified, direction 'down'
    #undef ALCH1
    #undef ALCH2

    #undef NORMAL
    #undef EXCLUDED
    #undef MODIFIED
    
    /********************************************************************/
    /******NONBONDEDBASE2 FOR EXCLUDED INTERACTIONS SCALED BY LAMBDA*****/
    /********************************************************************/
    #ifdef FULLELECT
    #undef FAST
    #define FAST(X)
    #define NORMAL(X)
    #define EXCLUDED(X) X
    #define MODIFIED(X)

    #define ALCH1(X) X
    #define ALCH2(X)
    npairi = pairlist_from_pairlist(ComputeNonbondedUtil::cutoff2,
            p_i_x, p_i_y, p_i_z, p_1, pairlistxA1_save, npairxA1, pairlisti,
            r2_delta, r2list);
        exclChecksum += npairi;
        SELF(
        for (k=0; k<npairi && pairlisti[k] < j_hgroup; ++k) --exclChecksum;
        )
    #include  "ComputeNonbondedBase2.h"  //excluded, direction 'up'
    #undef ALCH1
    #undef ALCH2

    #define ALCH1(X)
    #define ALCH2(X) X
    npairi = pairlist_from_pairlist(ComputeNonbondedUtil::cutoff2,
            p_i_x, p_i_y, p_i_z, p_1, pairlistxA2_save, npairxA2, pairlisti,
            r2_delta, r2list);
        exclChecksum += npairi;
        SELF(
        for (k=0; k<npairi && pairlisti[k] < j_hgroup; ++k) --exclChecksum;
        )
    #include  "ComputeNonbondedBase2.h"  //excluded, direction 'down'
    #undef ALCH1
    #undef ALCH2

    #undef FAST
    #ifdef SLOWONLY
      #define FAST(X)
    #else
      #define FAST(X) X
    #endif
    #undef NORMAL
    #undef EXCLUDED
    #undef MODIFIED
    #else
        exclChecksum += npairxA1 + npairxA2;
        SELF(
          for (k=0; k<npairxA1 && pairlistxA1_save[k] < j_hgroup; ++k) --exclChecksum;
          for (k=0; k<npairxA2 && pairlistxA2_save[k] < j_hgroup; ++k) --exclChecksum;
        )
    #endif

    #undef ALCHPAIR
    #define ALCHPAIR(X)
    #undef NOT_ALCHPAIR
    #define NOT_ALCHPAIR(X) X

#endif // end nonbondedbase2.h includes for alchemical pairlists

#ifdef NETWORK_PROGRESS
    CkNetworkProgress();
#endif

#ifdef ARCH_POWERPC
    //data cache block touch the position structure
    __dcbt ((void *) &(p_0[i+1]));
    __prefetch_by_load ((void *)&(groupCount));
#endif

    SHORT( FAST( f_0[i].x += f_i_x; ) )
    SHORT( FAST( f_0[i].y += f_i_y; ) )
    SHORT( FAST( f_0[i].z += f_i_z; ) )
    FULL( fullf_0[i].x += fullf_i_x; )
    FULL( fullf_0[i].y += fullf_i_y; )
    FULL( fullf_0[i].z += fullf_i_z; )

	// PAIR( iout << i << " " << i_upper << " end\n" << endi;)
  } // for i

  // PAIR(iout << "++++++++\n" << endi;)

#ifdef f_1
#undef f_1
#endif
#if ( SHORT( FAST( 1+ ) ) 0 )
#if ( SELF( 1+ ) 0 )
  {
    Force *patch_f_0 = params->ff[0];

#ifndef ARCH_POWERPC 
#pragma ivdep
#endif
    for ( int i = 0; i < i_upper; ++i ) {
      patch_f_0[i].x += f_0[i].x;
      patch_f_0[i].y += f_0[i].y;
      patch_f_0[i].z += f_0[i].z;
      virial_xx += f_0[i].x * p_0[i].position.x;
      virial_xy += f_0[i].x * p_0[i].position.y;
      virial_xz += f_0[i].x * p_0[i].position.z;
      virial_yy += f_0[i].y * p_0[i].position.y;
      virial_yz += f_0[i].y * p_0[i].position.z;
      virial_zz += f_0[i].z * p_0[i].position.z;
    }
  }
#endif
#endif
#ifdef fullf_1
#undef fullf_1
#endif
#if ( FULL( 1+ ) 0 )
#if ( SELF( 1+ ) 0 )
  {
    Force *patch_fullf_0 = params->fullf[0];
#ifndef ARCH_POWERPC 
#pragma ivdep
#endif
    for ( int i = 0; i < i_upper; ++i ) {
      patch_fullf_0[i].x += fullf_0[i].x;
      patch_fullf_0[i].y += fullf_0[i].y;
      patch_fullf_0[i].z += fullf_0[i].z;
      fullElectVirial_xx += fullf_0[i].x * p_0[i].position.x;
      fullElectVirial_xy += fullf_0[i].x * p_0[i].position.y;
      fullElectVirial_xz += fullf_0[i].x * p_0[i].position.z;
      fullElectVirial_yy += fullf_0[i].y * p_0[i].position.y;
      fullElectVirial_yz += fullf_0[i].y * p_0[i].position.z;
      fullElectVirial_zz += fullf_0[i].z * p_0[i].position.z;
    }
  }
#endif
#endif

  reduction[exclChecksumIndex] += exclChecksum;
  FAST
  (
  ENERGY( reduction[vdwEnergyIndex] += vdwEnergy; )
  SHORT
  (
  ENERGY( reduction[electEnergyIndex] += electEnergy; )
  )
  ALCH
  (
    TI(reduction[vdwEnergyIndex_ti_1] += vdwEnergy_ti_1;) 
    TI(reduction[vdwEnergyIndex_ti_2] += vdwEnergy_ti_2;) 
    FEP( reduction[vdwEnergyIndex_s] += vdwEnergy_s; )
  SHORT
  (
    FEP( reduction[electEnergyIndex_s] += electEnergy_s; )
    TI(reduction[electEnergyIndex_ti_1] += electEnergy_ti_1;) 
    TI(reduction[electEnergyIndex_ti_2] += electEnergy_ti_2;) 
  )
  )
  SHORT
  (
  reduction[virialIndex_XX] += virial_xx;
  reduction[virialIndex_XY] += virial_xy;
  reduction[virialIndex_XZ] += virial_xz;
  reduction[virialIndex_YX] += virial_xy;
  reduction[virialIndex_YY] += virial_yy;
  reduction[virialIndex_YZ] += virial_yz;
  reduction[virialIndex_ZX] += virial_xz;
  reduction[virialIndex_ZY] += virial_yz;
  reduction[virialIndex_ZZ] += virial_zz;
  )
  )
  FULL
  (
  ENERGY( reduction[fullElectEnergyIndex] += fullElectEnergy; )
  ALCH
  (
    FEP( reduction[fullElectEnergyIndex_s] += fullElectEnergy_s; )
    TI(reduction[fullElectEnergyIndex_ti_1] += fullElectEnergy_ti_1;) 
    TI(reduction[fullElectEnergyIndex_ti_2] += fullElectEnergy_ti_2;) 
  )
  reduction[fullElectVirialIndex_XX] += fullElectVirial_xx;
  reduction[fullElectVirialIndex_XY] += fullElectVirial_xy;
  reduction[fullElectVirialIndex_XZ] += fullElectVirial_xz;
  reduction[fullElectVirialIndex_YX] += fullElectVirial_xy;
  reduction[fullElectVirialIndex_YY] += fullElectVirial_yy;
  reduction[fullElectVirialIndex_YZ] += fullElectVirial_yz;
  reduction[fullElectVirialIndex_ZX] += fullElectVirial_xz;
  reduction[fullElectVirialIndex_ZY] += fullElectVirial_yz;
  reduction[fullElectVirialIndex_ZZ] += fullElectVirial_zz;
  )

  delete [] excl_flags_buff;

}

