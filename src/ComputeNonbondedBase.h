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

  INT(
  BigReal *pressureProfileReduction = params->pressureProfileReduction;
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

  register const BigReal cutoff2 = ( savePairlists ? plcutoff2 :
				ComputeNonbondedUtil:: cutoff2 );
  register const BigReal groupcutoff2 = ( savePairlists ? groupplcutoff2 :
				ComputeNonbondedUtil:: groupcutoff2 );
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
  const BigReal* const table_four = ComputeNonbondedUtil:: table_noshort;
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

  plint grouplist_std[1005];
  plint fixglist_std[1005];  // list of non-fixed groups if fixedAtomsOn
  plint goodglist_std[1005];
  plint pairlistx_std[1005];
  plint pairlistm_std[1005];
  plint pairlist_std[1005];
  plint pairlist2_std[1005];
  plint pairlisti_std[1005];

  plint *grouplist = grouplist_std;
  plint *fixglist = fixglist_std;
  plint *goodglist = goodglist_std;
  plint *pairlistx = pairlistx_std;
  plint *pairlistm = pairlistm_std;
  plint *pairlist = pairlist_std;
  plint *pairlist2 = pairlist2_std;
  plint *pairlistn_save;  int npairn;
  plint *pairlistx_save;  int npairx;
  plint *pairlistm_save;  int npairm;
  plint *pairlisti = ( j_upper >= 1000 ? new plint[j_upper+5] : pairlisti_std );

  int fixg_upper = 0;
  int g_upper = 0;

  if ( savePairlists || ! usePairlists ) {

  if ( j_upper >= 1000 ) {
    grouplist = new plint[j_upper+5];
    fixglist = new plint[j_upper+5];
    goodglist = new plint[j_upper+5];
    pairlistx = new plint[j_upper+5];
    pairlistm = new plint[j_upper+5];
    pairlist = new plint[j_upper+5];
    pairlist2 = new plint[j_upper+5];
  }

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

  SELF ( int64 pairCount = ( (i_upper-1) * (int64)j_upper ) / 2; )
  PAIR ( int64 pairCount = i_upper * (int64)j_upper; )
  int64 minPairCount = ( pairCount * params->minPart ) / params->numParts;
  int64 maxPairCount = ( pairCount * params->maxPart ) / params->numParts;
  pairCount = 0;

  for ( i = 0; i < (i_upper SELF(- 1)); ++i )
  {
	// PAIR( iout << i << " " << i_upper << " start\n" << endi;)
    const CompAtom &p_i = p_0[i];
    register const BigReal p_i_x = p_i.position.x;
    register const BigReal p_i_y = p_i.position.y;
    register const BigReal p_i_z = p_i.position.z;

    if ( p_i.hydrogenGroupSize ) {
      int64 opc = pairCount;
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

  if ( savePairlists || ! usePairlists ) {

    if ( ! savePairlists ) pairlists.reset();  // limit space usage

    const ExclusionCheck *exclcheck = mol->get_excl_check_for_atom(p_i.id);
    const int excl_min = exclcheck->min;
    const int excl_max = exclcheck->max;
    const char * const excl_flags = exclcheck->flags - excl_min;

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
    register plint *pli = pairlist + pairlistindex;

    {
      register plint *gli = goodglist;
      const plint *glist = ( groupfixed ? fixglist : grouplist );
      SELF( const int gl = ( groupfixed ? fixg_lower : g_lower ); )
      const int gu = ( groupfixed ? fixg_upper : g_upper );
      register int g = PAIR(0) SELF(gl);
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
    plint *pairlistn = pairlists.newlist(j_upper + 5 + 1 + 5) SELF(+ 1);
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
	if ( ( ! (atomfixed && p_1[j].atomFixed) ) && (r2 <= cutoff2) ) {
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
	if ( (! p_1[j].atomFixed) && (r2 <= cutoff2) ) {
          int atom2 = p_1[j].id;
          if ( atom2 >= excl_min && atom2 <= excl_max ) *(pli++) = j;
          else *(plin++) = j;
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
       int atom2 = p_1[j2].id;
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
	if (r2 <= cutoff2) {
          if ( atom2 >= excl_min && atom2 <= excl_max ) *(pli++) = j;
          else *(plin++) = j;
        }
        atom2 = p_1[j2].id;
       }
      }
    }
    int npair2 = pli - pairlist2;
    if ( npair2 ) pairlist2[npair2] = pairlist2[npair2-1];

    plint *plix = pairlistx;
    plint *plim = pairlistm;
    plint *pln = pairlistn;
    int k=0;
    SELF(
    for (; pln < plin && *pln < j_hgroup; ++pln) {
      *(plix++) = *pln;  --exclChecksum;
    }
    pairlistn_skip = pln - pairlistn;
    for (; k < npair2 && pairlist2[k] < j_hgroup; ++k) {
      *(plix++) = pairlist2[k];  --exclChecksum;
    }
    )
    for (; k < npair2; ++k ) {
      int j = pairlist2[k];
      int atom2 = p_1[j].id;
      int excl_flag = excl_flags[atom2];
      switch ( excl_flag ) {
      case 0:  *(plin++) = j;  break;
      case 1:  *(plix++) = j;  break;
      case 2:  *(plim++) = j;  break;
      }
    }
    exclChecksum += (plix - pairlistx);
    exclChecksum += (plim - pairlistm);

    npairn = plin - pln;
    pairlistn_save = pln;
    pairlistn_save[npairn] = npairn ? pairlistn_save[npairn-1] : -1;
    pairlists.newsize(plin - pairlistn SELF(+ 1) + 1);

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

	// PAIR( iout << i << " " << i_upper << " save\n" << endi;)
  } else { // if ( savePairlists || ! usePairlists )
	// PAIR( iout << i << " " << i_upper << " use\n" << endi;)

    pairlists.nextlist(&pairlistn_save,&npairn);  --npairn;
    //if ( npairn > 1000 )
//	iout << i << " " << i_upper << " " << npairn << " n\n" << endi;
    SELF(
    int pairlistn_skip = *pairlistn_save;
    pairlistn_save += (pairlistn_skip + 1);
    npairn -= (pairlistn_skip + 1);
    )
    pairlists.nextlist(&pairlistx_save,&npairx);  --npairx;
    //if ( npairx > 1000 )
//	iout << i << " " << i_upper << " " << npairx << " x\n" << endi;
    // exclChecksum += npairx;
    pairlists.nextlist(&pairlistm_save,&npairm);  --npairm;
    //if ( npairm > 1000 )
//	iout << i << " " << i_upper << " " << npairm << " m\n" << endi;
    // exclChecksum += npairm;

  } // if ( savePairlists || ! usePairlists )

    FEP( BigReal *lambda_table_i = lambda_table + 6 * p_i.partition; )

    LES( BigReal *lambda_table_i =
			lambda_table + (lesFactor+1) * p_i.partition; )


    const BigReal kq_i = COLOUMB * p_i.charge * scaling * dielectric_1;
    const LJTable::TableEntry * const lj_row =
		ljTable->table_row(mol->atomvdwtype(p_i.id));

    SHORT( FAST( BigReal & f_i_x = f_0[i].x; ) )
    SHORT( FAST( BigReal & f_i_y = f_0[i].y; ) )
    SHORT( FAST( BigReal & f_i_z = f_0[i].z; ) )
    FULL( BigReal & fullf_i_x = fullf_0[i].x; )
    FULL( BigReal & fullf_i_y = fullf_0[i].y; )
    FULL( BigReal & fullf_i_z = fullf_0[i].z; )

    int npairi;
    int k;

    npairi = pairlist_from_pairlist(ComputeNonbondedUtil:: cutoff2,
	p_i_x, p_i_y, p_i_z, p_1, pairlistn_save, npairn, pairlisti);

#define NORMAL(X) X
#define EXCLUDED(X)
#define MODIFIED(X)
#include  "ComputeNonbondedBase2.h"
#undef NORMAL
#undef EXCLUDED
#undef MODIFIED

    npairi = pairlist_from_pairlist(ComputeNonbondedUtil:: cutoff2,
	p_i_x, p_i_y, p_i_z, p_1, pairlistm_save, npairm, pairlisti);
    exclChecksum -= npairm - npairi;

#define NORMAL(X)
#define EXCLUDED(X)
#define MODIFIED(X) X
#include  "ComputeNonbondedBase2.h"
#undef NORMAL
#undef EXCLUDED
#undef MODIFIED

#ifdef FULLELECT
    npairi = pairlist_from_pairlist(ComputeNonbondedUtil:: cutoff2,
	p_i_x, p_i_y, p_i_z, p_1, pairlistx_save, npairx, pairlisti);
    exclChecksum -= npairx - npairi;

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
#endif

	// PAIR( iout << i << " " << i_upper << " end\n" << endi;)
  } // for i

  // PAIR(iout << "++++++++\n" << endi;)

  if (grouplist != grouplist_std) delete [] grouplist;
  if (fixglist != fixglist_std) delete [] fixglist;
  if (goodglist != goodglist_std) delete [] goodglist;
  if (pairlist != pairlist_std) delete [] pairlist;
  if (pairlist2 != pairlist2_std) delete [] pairlist2;
  if (pairlistx != pairlistx_std) delete [] pairlistx;
  if (pairlistm != pairlistm_std) delete [] pairlistm;
  if (pairlisti != pairlisti_std) delete [] pairlisti;

  reduction[exclChecksumIndex] += exclChecksum;
  FAST
  (
  ENERGY( reduction[vdwEnergyIndex] += vdwEnergy; )
  SHORT
  (
  ENERGY( reduction[electEnergyIndex] += electEnergy; )
  )
  FEP
  (
  ENERGY( reduction[vdwEnergyIndex_s] += vdwEnergy_s; )
  SHORT
  (
  ENERGY( reduction[electEnergyIndex_s] += electEnergy_s; )
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
  FEP
  (
  ENERGY( reduction[fullElectEnergyIndex_s] += fullElectEnergy_s; )
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

