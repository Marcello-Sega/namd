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
  #define SLOWONLYNAME(X) FULLELECTNAME( X ## _slow )
#else
  #define FAST(X) X
  #define SLOWONLYNAME(X) FULLELECTNAME( X )
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
#ifdef FEPFLAG
  #define FEPNAME(X) LAST( X ## _fep )
  #define FEP(X) X
#else
  #define FEPNAME(X) LAST( X )
  #define FEP(X)
#endif

#define LAST(X) X

// see if things are really messed up
SELF( PAIR( foo bar ) )

// ************************************************************
// function header
void ComputeNonbondedUtil :: NAME
  ( nonbonded *params )

// function body
{
  if ( ComputeNonbondedUtil::commOnly ) return;

  // speedup variables
  BigReal *reduction = params->reduction;

  // local variables
  int exclChecksum = 0;
  FAST
  (
  BigReal vdwEnergy = 0;
  BigReal electEnergy = 0;

  FEP
  (
  BigReal vdwEnergy_s = 0;
  BigReal electEnergy_s = 0;
  )
  
  BigReal virial_xx = 0;
  BigReal virial_xy = 0;
  BigReal virial_xz = 0;
  BigReal virial_yy = 0;
  BigReal virial_yz = 0;
  BigReal virial_zz = 0;
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
  const BigReal* const fast_table = ComputeNonbondedUtil:: fast_table;
  )
  FULL
  (
  const BigReal* const scor_table = ComputeNonbondedUtil:: scor_table;
  const BigReal* const slow_table = ComputeNonbondedUtil:: slow_table;
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
  const BigReal r2_delta_1 = ComputeNonbondedUtil:: r2_delta_1;

  const int i_upper = params->numAtoms[0];
  register const int j_upper = params->numAtoms[1];
  register int j;
  register int i;
  const CompAtom *p_0 = params->p[0];
  const CompAtom *p_1 = params->p[1];

  // check for all fixed atoms
  if ( fixedAtomsOn ) {
    register int all_fixed = 1;
    for ( i = 0; all_fixed && i < i_upper; ++i)
      all_fixed = p_0[i].atomFixed;
    PAIR
    (
    for ( j = 0; all_fixed && j < j_upper; ++j)
      all_fixed = p_1[j].atomFixed;
    )
    if ( all_fixed ) return;
  }

  SELF( int j_hgroup = 0; )
  int pairlistindex=0;
  int pairlistoffset=0;
  int pairlist_std[1005];  // pad 1 + 4 for nonbonded group runover
  int pairlist2_std[1005];  // pad 1 + 4 for nonbonded group runover
  int *const pairlist = (j_upper < 1000 ? pairlist_std : new int[j_upper+5]);
  int *const pairlist2 = (j_upper < 1000 ? pairlist2_std : new int[j_upper+5]);

  FAST
  (
    Force *f_0 = params->ff[0];
    Force *f_1 = params->ff[1];
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

    FAST( Force & f_i = f_0[i]; )
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
    const int groupfixed = ( p_i.groupFixed );

    // If patch divisions are not made by hydrogen groups, then
    // hydrogenGroupSize is set to 1 for all atoms.  Thus we can
    // carry on as if we did have groups - only less efficiently.
    // An optimization in this case is to not rebuild the temporary
    // pairlist but to include every atom in it.  This should be a
    // a very minor expense.

    register const CompAtom *p_j = p_1;
    SELF( p_j += i+1; )

    PAIR( j = 0; )
    SELF
    (
      if ( p_i.hydrogenGroupSize ) {
        // exclude child hydrogens of i
        j_hgroup = i + p_i.hydrogenGroupSize;
      }
      // add all child or sister hydrogens of i
      for ( j = i + 1; j < j_hgroup; ++j ) {
	pairlist[pairlistindex++] = j;
	p_j++;
      }
    )

    // add remaining atoms to pairlist via hydrogen groups
    register int *pli = pairlist + pairlistindex;

    if ( groupfixed ) { // tuned assuming most atoms fixed
      while ( j < j_upper )
	{
	register int hgs = ( p_j->nonbondedGroupIsAtom ? 1 :
                                                p_j->hydrogenGroupSize );
	if ( ! (p_j->groupFixed) )
	{
	  // use a slightly large cutoff to include hydrogens
	  if ( square(p_j->position.x-p_i_x,p_j->position.y-p_i_y,p_j->position.z-p_i_z) <= groupcutoff2 ) {
	    pli[0] = j;   // copy over the next four in any case
	    pli[1] = j+1;
	    pli[2] = j+2;
	    pli[3] = j+3; // assume hgs <= 4
	    pli += hgs;
	  }
	}
	j += hgs;
	p_j += hgs;
	} // for j
    } else SELF( if ( j < j_upper ) ) { // tuned assuming no fixed atoms
      register BigReal p_j_x = p_j->position.x;
      register BigReal p_j_y = p_j->position.y;
      register BigReal p_j_z = p_j->position.z;
      while ( j < j_upper ) {
	register int hgs = ( p_j->nonbondedGroupIsAtom ? 1 :
                                                p_j->hydrogenGroupSize );
	p_j += ( ( j + hgs < j_upper ) ? hgs : 0 );
	register BigReal r2 = p_i_x - p_j_x;
	r2 *= r2;
	p_j_x = p_j->position.x;	// preload
	register BigReal t2 = p_i_y - p_j_y;
	r2 += t2 * t2;
	p_j_y = p_j->position.y;	// preload
	t2 = p_i_z - p_j_z;
	r2 += t2 * t2;
	p_j_z = p_j->position.z;	// preload
	// use a slightly large cutoff to include hydrogens
	if ( r2 <= groupcutoff2 ) {
	  pli[0] = j;   // copy over the next four in any case
	  pli[1] = j+1;
	  pli[2] = j+2;
	  pli[3] = j+3; // assume hgs <= 4
	  pli += hgs;
	}
	j += hgs;
      } // for j
    }

    pairlistindex = pli - pairlist;
    // make sure padded element on pairlist points to real data
    if ( pairlistindex ) pairlist[pairlistindex] = pairlist[pairlistindex-1];
  } // if i is hydrogen group parent
  SELF
    (
    // self-comparisions require list to be incremented
    // pair-comparisions use entire list (pairlistoffset is 0)
    else pairlistoffset++;
    )

    const int atomfixed = ( p_i.atomFixed );

  FEP
    (
    const int ifep_type = p_i.partition;
    )

    const BigReal kq_i = COLOUMB * p_i.charge * scaling * dielectric_1;
    const LJTable::TableEntry * const lj_row =
		ljTable->table_row(mol->atomvdwtype(p_i.id));

    register int *pli = pairlist2;
    if ( atomfixed ) {
      for (int k=pairlistoffset; k<pairlistindex; k++) {
        j = pairlist[k];
        register const CompAtom *p_j = p_1 + j;
	if ( ! (p_j->atomFixed) ) {
          register BigReal r2 = square(p_j->position.x-p_i_x,
			p_j->position.y-p_i_y,p_j->position.z-p_i_z);
          if ( (r2 <= cutoff2) && ! ((r2 == 0) && ++exclChecksum) ) {
            *(pli++) = j;
          }
        }
      }
    } else {
      register const CompAtom *p_j = p_1;
      if ( pairlistoffset < pairlistindex ) p_j += pairlist[pairlistoffset];
      register BigReal p_j_x = p_j->position.x;
      register BigReal p_j_y = p_j->position.y;
      register BigReal p_j_z = p_j->position.z;
      int j1 = pairlist[pairlistoffset];

      for (int k=pairlistoffset; k<pairlistindex; k++) {
        int j = j1;

        // don't worry about [k+1] going beyond array since array is 1 too large
        j1 = pairlist[k+1];
        p_j = p_1 + j1;			// preload

	register BigReal r2 = p_i_x - p_j_x;
	r2 *= r2;
	p_j_x = p_j->position.x;	// preload
	register BigReal t2 = p_i_y - p_j_y;
	r2 += t2 * t2;
	p_j_y = p_j->position.y;	// preload
	t2 = p_i_z - p_j_z;
	r2 += t2 * t2;
	p_j_z = p_j->position.z;	// preload

        *pli = j;
        if ( (r2 <= cutoff2) && ! ((r2 == 0) && ++exclChecksum) ) { ++pli; }
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

      int table_i = (int) ( r2_delta_1 * r2 );
      FAST(
      const BigReal* const fast_i = fast_table + 4*table_i;
      BigReal fast_a = fast_i[0];
      )
      FULL(
      const BigReal* const scor_i = scor_table + 4*table_i;
      BigReal slow_a = scor_i[0]; 
      )

      const LJTable::TableEntry * lj_pars = 
		lj_row + 2 * mol->atomvdwtype(p_j->id);

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
      BigReal diffa = r2 - r2_delta * table_i;

    FEP
      (
      int jfep_type = p_j->partition;
      BigReal lambda_pair = 1.0;
      BigReal d_lambda_pair = 1.0;
      if (ifep_type || jfep_type) {
        if (ifep_type && jfep_type && ifep_type != jfep_type) {
	  lambda_pair = 0.0;
	  d_lambda_pair = 0.0;
        } else if ( fepOn ) {
          if (ifep_type == 1 || jfep_type == 1) {
	    lambda_pair = lambda;
	    d_lambda_pair = lambda2;
          } else if ( ifep_type == 2 || jfep_type == 2) {
	    lambda_pair = 1.0 - lambda;
	    d_lambda_pair = 1.0 - lambda2;
          }
        } else { // lesOn
          int fep_type = ifep_type ? ifep_type : jfep_type;
          lambda_pair = fep_type > lesFactor ? 0.0 : lesScaling;
	  d_lambda_pair = 0.0;
        }
      }
      )

      FAST
      (
      const BigReal A = scaling * lj_pars->A;
      const BigReal B = scaling * lj_pars->B;

      const BigReal r_2 = 1.0 / r2;
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

      vdwEnergy += FEP(lambda_pair *) switchVal * AmBterm;
      BigReal force_r = FEP(lambda_pair *)
        ( switchVal * 3.0 * (A*r_12 + AmBterm) * r_2 - AmBterm * dSwitchVal );

      FEP( vdwEnergy_s += d_lambda_pair * switchVal * AmBterm; )

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

      electEnergy += FEP(lambda_pair *) kqq * fast_val;
      force_r -= FEP(lambda_pair *) kqq * fast_dir;

      FEP( electEnergy_s += d_lambda_pair * kqq * fast_val; )

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

      fullElectEnergy += FEP(lambda_pair *) kqq * slow_val;
      BigReal fullforce_r = -1.0 * kqq * slow_dir FEP(* lambda_pair);
      
      FEP( fullElectEnergy_s += d_lambda_pair * kqq * slow_val; )

      {
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
  if (pairlist != pairlist_std) delete [] pairlist;
  if (pairlist2 != pairlist2_std) delete [] pairlist2;

  reduction[exclChecksumIndex] += exclChecksum;
  FAST
  (
  reduction[vdwEnergyIndex] += vdwEnergy;
  reduction[electEnergyIndex] += electEnergy;
  FEP
  (
  reduction[vdwEnergyIndex_s] += vdwEnergy_s;
  reduction[electEnergyIndex_s] += electEnergy_s;
  )
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

