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
#undef PPROF
#undef LAM
#define FEPNAME(X) LAST( X )
#define FEP(X)
#define LES(X)
#define INT(X)
#define PPROF(X)
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
  // const int r2_delta_expc = 64 * (r2_delta_exp - 127);
  const int r2_delta_expc = 64 * (r2_delta_exp - 1023);

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

  NBWORKARRAYSINIT(params->workArrays);

  int arraysize = j_upper+5;

  NBWORKARRAY(plint,pairlisti,arraysize)
  NBWORKARRAY(BigReal,r2list,arraysize)

  union { double f; int32 i[2]; } byte_order_test;
  byte_order_test.f = 1.0;  // should occupy high-order bits only
  int32 *r2iilist = (int32*)r2list + ( byte_order_test.i[0] ? 0 : 1 );

  if ( ! ( savePairlists || ! usePairlists ) ) arraysize = 0;

  NBWORKARRAY(plint,grouplist,arraysize);
  NBWORKARRAY(plint,fixglist,arraysize);
  NBWORKARRAY(plint,goodglist,arraysize);
  NBWORKARRAY(plint,pairlistx,arraysize);
  NBWORKARRAY(plint,pairlistm,arraysize);
  NBWORKARRAY(plint,pairlist,arraysize);
  NBWORKARRAY(plint,pairlist2,arraysize);

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
#if defined(ARCH_POWERPC) & !defined(MEM_OPT_VERSION)
	__dcbt((void *) &(p_0[i+1]));
#endif
        continue;
      }
    }

    register const BigReal p_i_x = p_i.position.x;
    register const BigReal p_i_y = p_i.position.y;
    register const BigReal p_i_z = p_i.position.z;

    PPROF(
        const int p_i_partition = p_i.partition;
        int n1 = (int)floor((p_i_z-pressureProfileMin)*invThickness);
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
      register plint *gli = goodglist;
      const plint *glist = ( groupfixed ? fixglist : grouplist );
      SELF( const int gl = ( groupfixed ? fixg_lower : g_lower ); )
      const int gu = ( groupfixed ? fixg_upper : g_upper );
      register int g = PAIR(0) SELF(gl);

      if ( g < gu ) {
	int hu = 0;
	if ( gu - g  >  6 ) { 

	  register  int jprev0 = glist[g];
	  register  int jprev1 = glist[g + 1];
	  
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
	    r2_0   =  t_0 * t_0 + r2_delta;
	    r2_1   =  t_1 * t_1 + r2_delta;
	    
	    t_0    =  p_i_y - pj_y_0;
	    t_1    =  p_i_y - pj_y_1;
	    r2_0  +=  t_0 * t_0;
	    r2_1  +=  t_1 * t_1;
	    
	    t_0    =  p_i_z - pj_z_0;
	    t_1    =  p_i_z - pj_z_1;
	    r2_0  +=  t_0 * t_0;
	    r2_1  +=  t_1 * t_1;
	    
	    jprev0     =  glist[g];
	    jprev1     =  glist[g+1];
	    
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
	    goodglist [ hu ] = j0;
	    goodglist [ hu + test0 ] = j1;
	    
	    hu += test0 + test1;
	  }
	  g-=2;
	}
	
	for (; g < gu; g++) {
	  int j = glist[g];
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
	    r2_0   =  t_0 * t_0 + r2_delta;
	    r2_1   =  t_1 * t_1 + r2_delta;
	    
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
    int k=0;
    SELF(
    for (; pln < plin && *pln < j_hgroup; ++pln) {
      *(plix++) = *pln;  // --exclChecksum;
    }
    pairlistn_skip = pln - pairlistn;
    for (; k < npair2 && pairlist2[k] < j_hgroup; ++k) {
      *(plix++) = pairlist2[k];  // --exclChecksum;
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
    // exclChecksum += (plix - pairlistx);
    // exclChecksum += (plim - pairlistm);

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

#ifdef NETWORK_PROGRESS
    CkNetworkProgress();
#endif

#if defined(ARCH_POWERPC) & !defined(MEM_OPT_VERSION)
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

  delete [] excl_flags_buff;

}

