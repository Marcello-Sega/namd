/**
***  Copyright (c) 1995, 1996, 1997, 1998, 1999, 2000 by
***  The Board of Trustees of the University of Illinois.
***  All rights reserved.
**/

// Several special cases are defined:
//   NBTYPE: exclusion method (NBPAIR, NBSELF -- mutually exclusive)
//   FULLELECT full electrostatics calculation?

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
  #define CLASSNAME(X) SLOWONLYNAME( X ## _pair )
#else
  #define PAIR(X)
#endif

#undef SELF
#if NBTYPE == NBSELF
  #define SELF(X) X
  #define CLASS ComputeNonbondedSelf
  #define CLASSNAME(X) SLOWONLYNAME( X ## _self )
#else
  #define SELF(X)
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

#undef FEPNAME
#undef FEP
#undef LES
#undef INT
#undef LAM
#define FEPNAME(X) LAST( X )
#define FEP(X)
#define LES(X)
#define INT(X)
#define LAM(X)
#ifdef FEPFLAG
  #undef FEPNAME
  #undef FEP
  #undef LAM
  #define FEPNAME(X) LAST( X ## _fep )
  #define FEP(X) X
  #define LAM(X) X
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

#define LAST(X) X

// see if things are really messed up
SELF( PAIR( foo bar ) )
LES( FEP( foo bar ) )
LES( INT( foo bar ) )
FEP( INT( foo bar ) )
LAM( INT( foo bar ) )

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

  // local variables
  int exclChecksum = 0;
  FAST
  (
  BigReal vdwEnergy = 0;
  SHORT
  (
  BigReal electEnergy = 0;
  )

  FEP
  (
  BigReal vdwEnergy_s = 0;
  SHORT
  (
  BigReal electEnergy_s = 0;
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
  BigReal fullElectEnergy = 0;
  FEP
  (
  BigReal fullElectEnergy_s = 0;
  )
  BigReal fullElectVirial_xx = 0;
  BigReal fullElectVirial_xy = 0;
  BigReal fullElectVirial_xz = 0;
  BigReal fullElectVirial_yy = 0;
  BigReal fullElectVirial_yz = 0;
  BigReal fullElectVirial_zz = 0;
  )

  // Bringing stuff into local namespace for speed.

  register const BigReal cutoff2 = ComputeNonbondedUtil:: cutoff2;
  register const BigReal groupcutoff2 = ComputeNonbondedUtil:: groupcutoff2;
  const BigReal dielectric_1 = ComputeNonbondedUtil:: dielectric_1;
  const LJTable* const ljTable = ComputeNonbondedUtil:: ljTable;
  LJTable::TableEntry ljNull;  ljNull.A = 0; ljNull.B = 0;
  const LJTable::TableEntry* const lj_null_pars = &ljNull;
  const Molecule* const mol = ComputeNonbondedUtil:: mol;
  FAST
  (
  SHORT
  (
  const BigReal* const fast_table = ComputeNonbondedUtil:: fast_table;
  )
  )
  FULL
  (
  SHORT
  (
  const BigReal* const scor_table = ComputeNonbondedUtil:: scor_table;
  const BigReal* const slow_table = ComputeNonbondedUtil:: slow_table;
  )
  NOSHORT
  (
  const BigReal* const scor_table = ComputeNonbondedUtil:: corr_table;
  const BigReal* const slow_table = ComputeNonbondedUtil:: full_table;
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
  const int r2_delta_expc = 64 * (r2_delta_exp - 127);

  const int i_upper = params->numAtoms[0];
  register const int j_upper = params->numAtoms[1];
  register int j;
  register int i;
  const CompAtom *p_0 = params->p[0];
  const CompAtom *p_1 = params->p[1];

  int grouplist_std[1005];
  int fixglist_std[1005];  // list of non-fixed groups if fixedAtomsOn
  int goodglist_std[1005];
  int *const grouplist = (j_upper < 1000 ? grouplist_std : new int[j_upper+5]);
  int *const fixglist = (j_upper < 1000 ? fixglist_std : new int[j_upper+5]);
  int *const goodglist = (j_upper < 1000 ? goodglist_std : new int[j_upper+5]);

  register int g = 0;
  for ( j = 0; j < j_upper; ++j ) {
    if ( p_1[j].hydrogenGroupSize || p_1[j].nonbondedGroupIsAtom ) {
      grouplist[g++] = j;
    }
  }
  const int g_upper = g;
  if ( g_upper ) grouplist[g_upper] = grouplist[g_upper-1];
  int fixg = 0;

  // check for all fixed atoms
  if ( fixedAtomsOn ) {
    register int all_fixed = 1;
    for ( g = 0; g < g_upper; ++g ) {
      j = grouplist[g];
      if ( ! p_1[j].groupFixed ) {
        all_fixed = 0;
        fixglist[fixg++] = j;
      }
    }
    PAIR
    (
    for ( i = 0; all_fixed && i < i_upper; ++i ) {
      if ( ! p_0[i].atomFixed ) all_fixed = 0;
    }
    )
    if ( all_fixed ) {
      if (grouplist != grouplist_std) delete [] grouplist;
      if (fixglist != fixglist_std) delete [] fixglist;
      if (goodglist != goodglist_std) delete [] goodglist;
      return;
    }
  }

  const int fixg_upper = fixg;
  if ( fixg_upper ) fixglist[fixg_upper] = fixglist[fixg_upper-1];

  SELF(
  int j_hgroup = 0;
  int g_lower = 0;
  int fixg_lower = 0;
  )
  int pairlistindex=0;
  int pairlistoffset=0;
  int pairlist_std[1005];  // pad 1 + 4 for nonbonded group runover
  int pairlist2_std[1005];  // pad 1 + 4 for nonbonded group runover
  int *const pairlist = (j_upper < 1000 ? pairlist_std : new int[j_upper+5]);
  int *const pairlist2 = (j_upper < 1000 ? pairlist2_std : new int[j_upper+5]);

  SHORT
  (
  FAST
  (
    Force *f_0 = params->ff[0];
    Force *f_1 = params->ff[1];
  )
  )
  FULL
  (
    Force *fullf_0 = params->fullf[0];
    Force *fullf_1 = params->fullf[1];
  )

  SELF ( int pairCount = ( (i_upper-1) * j_upper ) / 2; )
  PAIR ( int pairCount = i_upper * j_upper; )
  int minPairCount = ( pairCount * params->minPart ) / params->numParts;
  int maxPairCount = ( pairCount * params->maxPart ) / params->numParts;
  pairCount = 0;

  for ( i = 0; i < (i_upper SELF(- 1)); ++i )
  {
    const CompAtom &p_i = p_0[i];
    const ExclusionCheck *exclcheck = mol->get_excl_check_for_atom(p_i.id);
    const int excl_min = exclcheck->min;
    const int excl_max = exclcheck->max;
    const char * const excl_flags = exclcheck->flags;
    register const BigReal p_i_x = p_i.position.x;
    register const BigReal p_i_y = p_i.position.y;
    register const BigReal p_i_z = p_i.position.z;

    SHORT( FAST( Force & f_i = f_0[i]; ) )
    FULL( Force & fullf_i = fullf_0[i]; )

  if (p_i.hydrogenGroupSize || p_i.nonbondedGroupIsAtom)
    {
    if ( p_i.hydrogenGroupSize ) {
      int opc = pairCount;
      int hgs = p_i.hydrogenGroupSize;
      SELF
      (
      pairCount += hgs * ( i_upper - 1 - i );
      pairCount -= hgs * ( hgs - 1 ) / 2;
      )
      PAIR
      (
      pairCount += hgs * j_upper;
      )
      if ( opc < minPairCount || opc >= maxPairCount ) {
        i += hgs - 1;
        continue;
      }
    }

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
        j_hgroup = i + p_i.hydrogenGroupSize;
        while ( grouplist[g_lower] < j_hgroup ) ++g_lower;
        while ( fixglist[fixg_lower] < j_hgroup ) ++fixg_lower;
      }
      // add all child or sister hydrogens of i
      for ( j = i + 1; j < j_hgroup; ++j ) {
	pairlist[pairlistindex++] = j;
      }
    )

    // add remaining atoms to pairlist via hydrogen groups
    register int *pli = pairlist + pairlistindex;

    {
      register int *gli = goodglist;
      const int *glist = ( groupfixed ? fixglist : grouplist );
      SELF( const int gl = ( groupfixed ? fixg_lower : g_lower ); )
      const int gu = ( groupfixed ? fixg_upper : g_upper );
      g = PAIR(0) SELF(gl);
      if ( g < gu ) {
       int j2 = glist[g];
       BigReal p_j_x = p_1[j2].position.x;
       BigReal p_j_y = p_1[j2].position.y;
       BigReal p_j_z = p_1[j2].position.z;
       while ( g < gu ) {
        j = j2;
        j2 = glist[++g];
	BigReal r2 = p_i_x - p_j_x;
	r2 *= r2;
        p_j_x = p_1[j2].position.x;
	BigReal t2 = p_i_y - p_j_y;
	r2 += t2 * t2;
        p_j_y = p_1[j2].position.y;
	t2 = p_i_z - p_j_z;
	r2 += t2 * t2;
        p_j_z = p_1[j2].position.z;
	// use a slightly large cutoff to include hydrogens
	if ( r2 <= groupcutoff2 ) { *gli = j; ++gli; }
       }

       int hu = gli - goodglist;
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
    } PAIR( else {  // skip empty loops if no pairs were found
       int hgs = ( p_i.nonbondedGroupIsAtom ? 1 : p_i.hydrogenGroupSize );
       i += hgs - 1;
       continue;
    } )
  } // if i is hydrogen group parent
  SELF
    (
    // self-comparisions require list to be incremented
    // pair-comparisions use entire list (pairlistoffset is 0)
    else pairlistoffset++;
    )

    const int atomfixed = ( fixedAtomsOn && p_i.atomFixed );

    FEP( BigReal *lambda_table_i = lambda_table + 6 * p_i.partition; )

    LES( BigReal *lambda_table_i =
			lambda_table + (lesFactor+1) * p_i.partition; )


    const BigReal kq_i = COLOUMB * p_i.charge * scaling * dielectric_1;
    const LJTable::TableEntry * const lj_row =
		ljTable->table_row(mol->atomvdwtype(p_i.id));

    register int *pli = pairlist2;

    INT(
    if ( 1 ) {
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
	if ( ( ! (atomfixed && p_1[j].atomFixed) ) &&
	     (r2 <= cutoff2) && ! ((r2 <= r2_delta) && ++exclChecksum) ) {
          *(pli++) = j;
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
	if ( (! p_1[j].atomFixed) &&
	     (r2 <= cutoff2) && ! ((r2 <= r2_delta) && ++exclChecksum) ) {
          *(pli++) = j;
        }
      }
    } else {
      int k = pairlistoffset;
      int ku = pairlistindex;
      if ( k < ku ) {
       int j2 = pairlist[k];
       BigReal p_j_x = p_1[j2].position.x;
       BigReal p_j_y = p_1[j2].position.y;
       BigReal p_j_z = p_1[j2].position.z;
       while ( k < ku ) {
        j = j2;
        j2 = pairlist[++k];
	BigReal r2 = p_i_x - p_j_x;
	r2 *= r2;
        p_j_x = p_1[j2].position.x;
	BigReal t2 = p_i_y - p_j_y;
	r2 += t2 * t2;
        p_j_y = p_1[j2].position.y;
	t2 = p_i_z - p_j_z;
	r2 += t2 * t2;
        p_j_z = p_1[j2].position.z;
	if ( (r2 <= cutoff2) && ! ((r2 <= r2_delta) && ++exclChecksum) ) {
          *(pli++) = j;
        }
       }
      }
    }
    int npair2 = pli - pairlist2;

    for (int k=0; k<npair2; ++k) {

      const int j = pairlist2[k];
      register const CompAtom *p_j = p_1 + j;

      register const BigReal p_ij_x = p_i_x - p_j->position.x;
      register BigReal r2 = p_ij_x * p_ij_x;
      register const BigReal p_ij_y = p_i_y - p_j->position.y;
      r2 += p_ij_y * p_ij_y;
      register const BigReal p_ij_z = p_i_z - p_j->position.z;
      r2 += p_ij_z * p_ij_z;

      float r2f = r2;
      const int table_i = ((*((int32 *)&r2f)) >> 17) + r2_delta_expc;

      FAST(
      const BigReal r_2 = 1.0 / r2;
      )

      const LJTable::TableEntry * lj_pars = 
		lj_row + 2 * mol->atomvdwtype(p_j->id);

      FAST(
      SHORT(
      const BigReal* const fast_i = fast_table + 4*table_i;
      BigReal fast_a = fast_i[0];
      )
      )
      FULL(
      const BigReal* const scor_i = scor_table + 4*table_i;
      BigReal slow_a = scor_i[0]; 
      )

      *((int32 *)&r2f) &= 0xfffe0000;

      BigReal modf = 0.0;
      int atom2 = p_j->id;
      register char excl_flag = ( (atom2 >= excl_min && atom2 <= excl_max) ?
					excl_flags[atom2-excl_min] : 0 );
      if ( excl_flag ) { ++exclChecksum; }
      SELF( if ( j < j_hgroup ) { excl_flag = EXCHCK_FULL; } )
      if ( excl_flag ) {
	if ( excl_flag == EXCHCK_FULL ) {
	  lj_pars = lj_null_pars;
	  modf = 1.0;
	} else {
	  ++lj_pars;
	  modf = modf_mod;
	}
      }

      BigReal kqq = kq_i * p_j->charge;
      const BigReal diffa = r2 - r2f;

      FEP(
      int jfep_type = p_j->partition;
      BigReal lambda_pair = lambda_table_i[2*jfep_type];
      BigReal d_lambda_pair = lambda_table_i[2*jfep_type+1];
      )

      LES( BigReal lambda_pair = lambda_table_i[p_j->partition]; )

      FAST
      (
      const BigReal A = scaling * lj_pars->A;
      const BigReal B = scaling * lj_pars->B;

      const BigReal r_6 = r_2*r_2*r_2;
      const BigReal r_12 = r_6*r_6;

      // Lennard-Jones switching function
      const BigReal c2 = cutoff2-r2;
      const BigReal c4 = c2*(cutoff2+2.0*r2-3.0*switchOn2);
      const BigReal switchVal =		// used for Lennard-Jones
			( r2 > switchOn2 ? c2*c4*c1 : 1.0 );
      const BigReal dSwitchVal =	// d switchVal / d r2
			( r2 > switchOn2 ? 0.5*c3*(c2*c2-c4) : 0.0 );

      const BigReal AmBterm = (A*r_6 - B) * r_6;

      vdwEnergy += LAM(lambda_pair *) switchVal * AmBterm;
      BigReal force_r = LAM(lambda_pair *)
        ( switchVal * 3.0 * (A*r_12 + AmBterm) * r_2 - AmBterm * dSwitchVal );

      FEP( vdwEnergy_s += d_lambda_pair * switchVal * AmBterm; )
      
      INT( 
      reduction[pairVDWForceIndex_X] += 2.0 * force_r * p_ij_x;
      reduction[pairVDWForceIndex_Y] += 2.0 * force_r * p_ij_y;
      reduction[pairVDWForceIndex_Z] += 2.0 * force_r * p_ij_z;
      )

      SHORT(
      BigReal modfc = 1.0 - modf;
      fast_a *= modfc;
      BigReal fast_d = modfc * fast_i[3];
      BigReal fast_c = modfc * fast_i[2];
      BigReal fast_b = modfc * fast_i[1];

      {
      register BigReal fast_val =
	( ( diffa * fast_d + fast_c ) * diffa + fast_b ) * diffa + fast_a;
      register BigReal fast_dir =
	( 3.0 * diffa * fast_d + 2.0 * fast_c ) * diffa + fast_b;

      electEnergy += LAM(lambda_pair *) kqq * fast_val;
      force_r -= LAM(lambda_pair *) kqq * fast_dir;

      FEP( electEnergy_s += d_lambda_pair * kqq * fast_val; )

      INT(
      reduction[pairElectForceIndex_X] -= 2.0 * kqq * fast_dir * p_ij_x;
      reduction[pairElectForceIndex_Y] -= 2.0 * kqq * fast_dir * p_ij_y;
      reduction[pairElectForceIndex_Z] -= 2.0 * kqq * fast_dir * p_ij_z;
      )
/*
      // JCP ERROR CHECKING CODE
      if ( abs(fast_val) > 1.0e4 || abs(fast_dir) > 1.0e4 ) {
        iout << iERROR << "atom " << p_i.id << " pos " << p_i.position << "\n";
        iout << iERROR << "atom " << p_j->id << " pos " << p_j->position << "\n";
        iout << iERROR << "r2 " << r2 << " h " << r2_delta << " i " << table_i << " d " << diffa << "\n";
        iout << iERROR << "fast " << fast_a << " " << fast_b << " " << fast_c << " " << fast_d  << "\n";
        iout << iERROR << "val " << fast_val << " dir " << fast_dir << " kqq " << kqq << "\n";
        iout << endi;
      }
*/
      }

      force_r *= 2.0;
      Force & f_j = f_1[j];
      register BigReal tmp_x = force_r * p_ij_x;
      virial_xx += tmp_x * p_ij_x;
      virial_xy += tmp_x * p_ij_y;
      virial_xz += tmp_x * p_ij_z;
      f_i.x += tmp_x;
      f_j.x -= tmp_x;
      register BigReal tmp_y = force_r * p_ij_y;
      virial_yy += tmp_y * p_ij_y;
      virial_yz += tmp_y * p_ij_z;
      f_i.y += tmp_y;
      f_j.y -= tmp_y;
      register BigReal tmp_z = force_r * p_ij_z;
      virial_zz += tmp_z * p_ij_z;
      f_i.z += tmp_z;
      f_j.z -= tmp_z;

      )
      )

      FULL(
      const BigReal* const slow_i = ( modf ? slow_table + 4*table_i : scor_i );
      slow_a -= modf * slow_i[0];
      BigReal slow_b = scor_i[1] - modf * slow_i[1];
      BigReal slow_c = scor_i[2] - modf * slow_i[2];
      BigReal slow_d = scor_i[3] - modf * slow_i[3];

      register BigReal slow_val =
	( ( diffa * slow_d + slow_c ) * diffa + slow_b ) * diffa + slow_a;
      register BigReal slow_dir =
	( 3.0 * diffa * slow_d + 2.0 * slow_c ) * diffa + slow_b;

      fullElectEnergy += LAM(lambda_pair *) kqq * slow_val;
      BigReal fullforce_r = -1.0 * kqq * slow_dir LAM(* lambda_pair);
      
      FEP( fullElectEnergy_s += d_lambda_pair * kqq * slow_val; )

      INT(
      reduction[pairElectForceIndex_X] -= 2.0 * kqq * slow_dir * p_ij_x;
      reduction[pairElectForceIndex_Y] -= 2.0 * kqq * slow_dir * p_ij_y;
      reduction[pairElectForceIndex_Z] -= 2.0 * kqq * slow_dir * p_ij_z;
      )

      {
      NOSHORT(
      fullforce_r += force_r;
      )
      fullforce_r *= 2.0;
      Force & fullf_j = fullf_1[j];
      register BigReal tmp_x = fullforce_r * p_ij_x;
      fullElectVirial_xx += tmp_x * p_ij_x;
      fullElectVirial_xy += tmp_x * p_ij_y;
      fullElectVirial_xz += tmp_x * p_ij_z;
      fullf_i.x += tmp_x;
      fullf_j.x -= tmp_x;
      register BigReal tmp_y = fullforce_r * p_ij_y;
      fullElectVirial_yy += tmp_y * p_ij_y;
      fullElectVirial_yz += tmp_y * p_ij_z;
      fullf_i.y += tmp_y;
      fullf_j.y -= tmp_y;
      register BigReal tmp_z = fullforce_r * p_ij_z;
      fullElectVirial_zz += tmp_z * p_ij_z;
      fullf_i.z += tmp_z;
      fullf_j.z -= tmp_z;

      }
      )

    } // for pairlist
  } // for i
  if (grouplist != grouplist_std) delete [] grouplist;
  if (fixglist != fixglist_std) delete [] fixglist;
  if (goodglist != goodglist_std) delete [] goodglist;
  if (pairlist != pairlist_std) delete [] pairlist;
  if (pairlist2 != pairlist2_std) delete [] pairlist2;

  reduction[exclChecksumIndex] += exclChecksum;
  FAST
  (
  reduction[vdwEnergyIndex] += vdwEnergy;
  SHORT
  (
  reduction[electEnergyIndex] += electEnergy;
  )
  FEP
  (
  reduction[vdwEnergyIndex_s] += vdwEnergy_s;
  SHORT
  (
  reduction[electEnergyIndex_s] += electEnergy_s;
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
  reduction[fullElectEnergyIndex] += fullElectEnergy;
  FEP
  (
  reduction[fullElectEnergyIndex_s] += fullElectEnergy_s;
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
}

