/**
***  Copyright (c) 1995, 1996, 1997, 1998, 1999, 2000 by
***  The Board of Trustees of the University of Illinois.
***  All rights reserved.
**/

// Several special cases are defined:
//   NBTYPE: exclusion method (NBPAIR, NBSELF, NBEXCL -- mutually exclusive)
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

#undef EXCL
#undef NOEXCL
#if NBTYPE == NBEXCL
  #define EXCL(X) X
  #define CLASS ComputeNonbondedExcl
  #define CLASSNAME(X) SLOWONLYNAME( X ## _excl )
  #define NOEXCL(X)
#else
  #define EXCL(X)
  #define NOEXCL(X) X
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
  #if FULLELECT == FULLELECT_PME
    #define FULLELECTNAME(X) SPLITTINGNAME( X ## _fullelect_pme )
  #else
    #define FULLELECTNAME(X) SPLITTINGNAME( X ## _fullelect )
  #endif
  #define FULL(X) X
  #define NOFULL(X)
#else
  #define FULLELECTNAME(X) SPLITTINGNAME( X )
  #define FULL(X)
  #define NOFULL(X) X
#endif

#undef SPLITTINGNAME
#if SPLIT_TYPE == SPLIT_NONE
  #define SPLITTINGNAME(X) LAST( X )
#elif SPLIT_TYPE == SPLIT_XPLOR
  #define SPLITTINGNAME(X) LAST( X ## _xplor )
#elif SPLIT_TYPE == SPLIT_C1
  #define SPLITTINGNAME(X) LAST( X ## _c1 )
#endif

#define LAST(X) X

// see if things are really messed up
SELF( PAIR( foo bar ) )
SELF( EXCL( foo bar ) )
EXCL( PAIR( foo bar ) )

// ************************************************************
// function header
void ComputeNonbondedUtil :: NAME
  ( nonbonded *params )

// function body
{
  // speedup variables
  BigReal *reduction = params->reduction;
  EXCL
  (
    Position & p_ij = params->p_ij;
    FAST
    (
    Force & f_i = params->ff[0][0];
    Force & f_j = params->ff[1][0];
    )
    const AtomProperties & a_i = params->a[0][0];
    const AtomProperties & a_j = params->a[1][0];
    int & m14 = params->m14;
    FULL
    (
    Force & fullf_i = params->fullf[0][0];
    Force & fullf_j = params->fullf[1][0];
    )
  )

  // local variables
  BigReal vdwEnergy = 0;
  BigReal electEnergy = 0;
  FAST
  (
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
  BigReal fullElectVirial_xx = 0;
  BigReal fullElectVirial_xy = 0;
  BigReal fullElectVirial_xz = 0;
  BigReal fullElectVirial_yy = 0;
  BigReal fullElectVirial_yz = 0;
  BigReal fullElectVirial_zz = 0;
  )
  NOEXCL
  (
  // Bringing stuff into local namespace for speed.
  // Probably makes things slower in exclusion mode, though.

  register const BigReal cutoff2 = ComputeNonbondedUtil:: cutoff2;
  register const BigReal groupcutoff2 = ComputeNonbondedUtil:: groupcutoff2;
  const BigReal dielectric_1 = ComputeNonbondedUtil:: dielectric_1;

  const LJTable* const ljTable = ComputeNonbondedUtil:: ljTable;
  const Molecule* const mol = ComputeNonbondedUtil:: mol;
  const BigReal scaling = ComputeNonbondedUtil:: scaling;
  const BigReal scale14 = ComputeNonbondedUtil:: scale14;
  const Real switchOn = ComputeNonbondedUtil:: switchOn;
  const BigReal switchOn2 = ComputeNonbondedUtil:: switchOn2;
  // const BigReal c0 = ComputeNonbondedUtil:: c0;
  const BigReal c1 = ComputeNonbondedUtil:: c1;
  const BigReal c3 = ComputeNonbondedUtil:: c3;
  const BigReal c5 = ComputeNonbondedUtil:: c5;
  const BigReal c6 = ComputeNonbondedUtil:: c6;
  const BigReal c7 = ComputeNonbondedUtil:: c7;
  const BigReal c8 = ComputeNonbondedUtil:: c8;
  // const BigReal d0 = ComputeNonbondedUtil:: d0;
  )

  BigReal kqq;	// initialized later
  BigReal f;	// initialized later

EXCL
(
  const BigReal kq_i = COLOUMB * a_i.charge * scaling * dielectric_1;
  register const BigReal p_ij_x = p_ij.x;
  register const BigReal p_ij_y = p_ij.y;
  register const BigReal p_ij_z = p_ij.z;
)

NOEXCL
(
  const int i_upper = params->numAtoms[0];
  register const int j_upper = params->numAtoms[1];
  register int j;
  register int i;

  // check for all fixed atoms
  if ( fixedAtomsOn ) {
    AtomProperties *a_i = params->a[0];
    register int all_fixed = 1;
    for ( i = 0; all_fixed && i < i_upper; ++i, ++a_i)
      all_fixed = a_i->flags & ATOM_FIXED;
    PAIR
    (
    AtomProperties *a_j = params->a[1];
    for ( j = 0; all_fixed && j < j_upper; ++j, ++a_j)
      all_fixed = a_j->flags & ATOM_FIXED;
    )
    if ( all_fixed ) return;
  }

  SELF( int j_hgroup = 0; )
  int pairlistindex=0;
  static int pairlist_std[1001];
  int pairlistoffset=0;
  int *pairlist = pairlist_std;

  // generate mini-pairlist.
  if (1000 < j_upper)
	{
	if (pairlist != pairlist_std) delete [] pairlist;
	pairlist = new int[j_upper + 1];
	pairlist[j_upper] = 0;
	}
  else
	{
	pairlist = pairlist_std;
	pairlist[1000] = 0;
	}

  // for speeding up the for-loop
  const AtomProperties *a_0 = params->a[0];
  const AtomProperties *a_1 = params->a[1];
  const Position *p_0 = params->p[0];
  const Position *p_1 = params->p[1];
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

  SELF
  (
  int pairCount = ( (i_upper-1) * j_upper ) / 2;
  int minPairCount = ( pairCount * params->minPart ) / params->numParts;
  int maxPairCount = ( pairCount * params->maxPart ) / params->numParts;
  pairCount = 0;
  )

  for ( i = 0; i < (i_upper SELF(- 1)); ++i )
  {
    const AtomProperties &a_i = a_0[i];
    const Position &p_i = p_0[i];
    register const BigReal p_i_x = p_i.x;
    register const BigReal p_i_y = p_i.y;
    register const BigReal p_i_z = p_i.z;

    FAST( Force & f_i = f_0[i]; )
    FULL( Force & fullf_i = fullf_0[i]; )

  if (a_i.nonbondedGroupSize) // if hydrogen group parent
    {
    SELF
    (
    if ( a_i.hydrogenGroupSize ) {
      int opc = pairCount;
      int hgs = a_i.hydrogenGroupSize;
      pairCount += hgs * ( i_upper - 1 - i );
      pairCount -= hgs * ( hgs - 1 ) / 2;
      if ( opc < minPairCount || opc >= maxPairCount ) {
        i += hgs - 1;
        continue;
      }
    }
    )

    pairlistindex = 0;	// initialize with 0 elements
    pairlistoffset=0;
    const int groupfixed = ( a_i.flags & GROUP_FIXED );

    // If patch divisions are not made by hydrogen groups, then
    // nonbondedGroupSize is set to 1 for all atoms.  Thus we can
    // carry on as if we did have groups - only less efficiently.
    // An optimization in this case is to not rebuild the temporary
    // pairlist but to include every atom in it.  This should be a
    // a very minor expense.

    register const Position *p_j = p_1;
    SELF( p_j += i+1; )

    PAIR( j = 0; )
    SELF
    (
      if ( a_i.hydrogenGroupSize ) {
        // exclude child hydrogens of i
        j_hgroup = i + a_i.hydrogenGroupSize;
      }
      // add all child or sister hydrogens of i
      for ( j = i + 1; j < j_hgroup; ++j ) {
	pairlist[pairlistindex++] = j;
	p_j++;
      }
    )

    // add remaining atoms to pairlist via hydrogen groups
    register const AtomProperties *pa_j = a_1 + j;
    register int *pli = pairlist + pairlistindex;

    if ( groupfixed ) { // tuned assuming most atoms fixed
      while ( j < j_upper )
	{
	register int hgs = pa_j->nonbondedGroupSize;
	if ( ! (pa_j->flags & GROUP_FIXED) )
	{
	  p_j = p_1 + j;
	  // use a slightly large cutoff to include hydrogens
	  if ( square(p_j->x-p_i_x,p_j->y-p_i_y,p_j->z-p_i_z) <= groupcutoff2 )
		{
		register int l = j;
		register int lm = j + hgs;
		for( ; l<lm; ++l)
		  {
		  *pli = l;
		  ++pli;
		  }
		}
	}
	j += hgs;
	pa_j += hgs;
	} // for j
    } else SELF( if ( j < j_upper ) ) { // tuned assuming no fixed atoms
      register BigReal p_j_x = p_j->x;
      register BigReal p_j_y = p_j->y;
      register BigReal p_j_z = p_j->z;
      while ( j < j_upper )
	{
	register int hgs = pa_j->nonbondedGroupSize;
	p_j += ( ( j + hgs < j_upper ) ? hgs : 0 );
	register BigReal r2 = p_i_x - p_j_x;
	r2 *= r2;
	p_j_x = p_j->x;					// preload
	register BigReal t2 = p_i_y - p_j_y;
	r2 += t2 * t2;
	p_j_y = p_j->y;					// preload
	t2 = p_i_z - p_j_z;
	r2 += t2 * t2;
	p_j_z = p_j->z;					// preload
	// use a slightly large cutoff to include hydrogens
	if ( r2 <= groupcutoff2 )
		{
		register int l = j;
		j += hgs;
		for( ; l<j; ++l)
		  {
		  *pli = l;
		  ++pli;
		  }
		}
	else j += hgs;
	pa_j += hgs;
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

    const int atomfixed = ( a_i.flags & ATOM_FIXED );

    const BigReal kq_i_u = COLOUMB * a_i.charge * scaling * dielectric_1;
    const BigReal kq_i_s = kq_i_u * scale14;

    register const Position *p_j = p_1;

    if ( pairlistoffset < pairlistindex ) p_j += pairlist[pairlistoffset];

    register BigReal p_j_x = p_j->x;
    register BigReal p_j_y = p_j->y;
    register BigReal p_j_z = p_j->z;

    for (int k=pairlistoffset; k<pairlistindex; k++)
    {
      j = pairlist[k];
      // don't worry about [k+1] going beyond array since array is 1 too large
      p_j += pairlist[k+1]-j; // preload
      register const BigReal p_ij_x = p_i_x - p_j_x;
      p_j_x = p_j->x;					// preload
      register const BigReal p_ij_y = p_i_y - p_j_y;
      p_j_y = p_j->y;					// preload
      register const BigReal p_ij_z = p_i_z - p_j_z;
      p_j_z = p_j->z;					// preload

)

      // common code
      register BigReal r2 = square(p_ij_x,p_ij_y,p_ij_z);

      if ( r2 > cutoff2 )
      {
	NOEXCL( continue; )
	EXCL(
	  FULL(
	    // Do a quick fix and get out!
	    const BigReal r_1 = 1.0/sqrt(r2);
	    kqq = kq_i * a_j.charge;
	    f = kqq*r_1;
	    if ( m14 ) f *= ( 1.0 - scale14 );
	    fullElectEnergy -= f;
	    const Vector f_elec = p_ij * ( f * r_1 * r_1 );
	    fullf_i -= f_elec;
	    fullf_j += f_elec;
	    reduction[fullElectEnergyIndex] += fullElectEnergy;
	    reduction[fullElectVirialIndex_XX] -= f_elec.x * p_ij.x;
	    reduction[fullElectVirialIndex_XY] -= f_elec.x * p_ij.y;
	    reduction[fullElectVirialIndex_XZ] -= f_elec.x * p_ij.z;
	    reduction[fullElectVirialIndex_YX] -= f_elec.y * p_ij.x;
	    reduction[fullElectVirialIndex_YY] -= f_elec.y * p_ij.y;
	    reduction[fullElectVirialIndex_YZ] -= f_elec.y * p_ij.z;
	    reduction[fullElectVirialIndex_ZX] -= f_elec.z * p_ij.x;
	    reduction[fullElectVirialIndex_ZY] -= f_elec.z * p_ij.y;
	    reduction[fullElectVirialIndex_ZZ] -= f_elec.z * p_ij.z;
	  )
	return; )
      }

NOEXCL
(
      const AtomProperties & a_j = a_1[j];

      if ( atomfixed && ( a_j.flags & ATOM_FIXED ) ) continue;

      BigReal kq_i = kq_i_u;

      FULL
      (
        Force & fullf_j = fullf_1[j];
	const BigReal r = sqrt(r2);
        const BigReal r_1 = 1.0/r;
        kqq = kq_i * a_j.charge;
        f = kqq*r_1;

	//  Patch code for full electrostatics algorithms which
	//  need to do a direct component like PME or P3M fits here.
	//  This may not be the absolutely most efficient position,
	//  but it is nicely detached from the complex logic of the
	//  exclusion checking system below so it won't break anything.

	switch ( FULLELECT ) {  // compiler should optimize this away
	  case FULLELECT_PME: {
	    register BigReal tmp_a = r * ewaldcof;
	    register BigReal tmp_b = erfc(tmp_a);
	    register BigReal tmp_c =
		pi_ewaldcof*exp(-(tmp_a*tmp_a))*r_1 + tmp_b*r_1*r_1;
	    register BigReal tmp_f = tmp_c * f;
	    fullElectEnergy += tmp_b * f;
	    register BigReal tmp_x = tmp_f * p_ij_x;
	    fullElectVirial_xx += tmp_x * p_ij_x;
	    fullElectVirial_xy += tmp_x * p_ij_y;
	    fullElectVirial_xz += tmp_x * p_ij_z;
	    fullf_i.x += tmp_x;
	    fullf_j.x -= tmp_x;
	    register BigReal tmp_y = tmp_f * p_ij_y;
	    fullElectVirial_yy += tmp_y * p_ij_y;
	    fullElectVirial_yz += tmp_y * p_ij_z;
	    fullf_i.y += tmp_y;
	    fullf_j.y -= tmp_y;
	    register BigReal tmp_z = tmp_f * p_ij_z;
	    fullElectVirial_zz += tmp_z * p_ij_z;
	    fullf_i.z += tmp_z;
	    fullf_j.z -= tmp_z;
	  } break;
	}
      )
)

      register BigReal force_r = 0.0;		//  force / r
      FULL
      (
      register BigReal fullforce_r = 0.0;	//  fullforce / r
      )

      const LJTable::TableEntry * lj_pars = 
		ljTable->table_val(a_i.type, a_j.type);

      if ( r2 <= lj_pars->exclcut2 SELF( || j < j_hgroup ) )
      {
	NOEXCL
	(
	register int *list;
	register int other_int;
	SELF( if ( ! ( j < j_hgroup ) ) ) {
           //  We want to search the array of the smaller atom
	  int atom1 = a_i.id;
	  int atom2 = a_j.id;
	  list = mol->get_excl_check_for_atom(atom1<atom2?atom1:atom2);
	  other_int = atom1<atom2?atom2:atom1;
	  for ( ; *list != other_int && *list != -1; ++list );
	}
	if ( SELF( j < j_hgroup || ) ! ( *list != other_int ) ) {
	  FULL
	  (
	    // Do a quick fix and get out!
	    fullElectEnergy -= f;
	    fullforce_r = -f * r_1*r_1;
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
	  )
	  continue;  // Must have stored force by now.
	} else {
	  for ( ++list ; *list != other_int && *list != -1; ++list );
          if ( ! ( *list != other_int ) ) {
	    FULL
	    (
	      // Make full electrostatics match rescaled charges!
	      f *= ( 1. - scale14 );
	      fullElectEnergy -= f;
	      fullforce_r -= f * r_1*r_1;
	    )
	    lj_pars = ljTable->table_val(a_i.type, a_j.type, 1);
	    kq_i = kq_i_s;
	  }
	}
	)

        EXCL
        (
	  return;
        )

      }

      NOEXCL
      (
      FAST( Force & f_j = f_1[j]; )
      )

      NOFULL
      (
      const BigReal r = sqrt(r2);
      const BigReal r_1 = 1.0/r;
      )

      FULL
      (
        EXCL
        (
            const BigReal r = sqrt(r2);
            const BigReal r_1 = 1.0/r;
        )
      )

      BigReal switchVal; // used for Lennard-Jones
      BigReal shiftVal; // used for electrostatics splitting as well
      BigReal dSwitchVal; // used for Lennard-Jones
      BigReal dShiftVal; // used for electrostatics splitting as well

      // Lennard-Jones switching function
      if (r > switchOn)
      {
	const BigReal c2 = cutoff2-r2;
	const BigReal c4 = c2*(cutoff2+2.0*r2-3.0*switchOn2);
	switchVal = c2*c4*c1;
	dSwitchVal = c3*r*(c2*c2-c4);
      }
      else
      {
	switchVal = 1;
	dSwitchVal = 0;
      }


//  --------------------------------------------------------------------------
//  BEGIN SHIFTING / SPLITTING FUNCTION DEFINITIONS
//  --------------------------------------------------------------------------


    // compiler should optimize this code easily
    switch(SPLIT_TYPE)
      {
      case SPLIT_NONE:
	// shifting(shiftVal,dShiftVal,r,r2,c5,c6);
	shiftVal = 1.0 - r2*c5;
	dShiftVal = c6*shiftVal*r;
	shiftVal *= shiftVal;
	break;
      case SPLIT_C1:
	{
	  const BigReal d1 = r2 * c7;
	  shiftVal = 1.0 + r * ( d1 - c8 );
	  dShiftVal = 3.0 * d1 - c8;
	}
	/*
	// c1splitting(shiftVal,dShiftVal,r,d0,switchOn);
	dShiftVal = 0;  // formula only correct for forces
	if (r > switchOn) {
	  const BigReal d1 = d0*(r-switchOn);
	  shiftVal = 1.0 + d1*d1*(2.0*d1-3.0);
	} else {
	  shiftVal = 1;
	}
	*/
	break;
      case SPLIT_XPLOR:
	// xplorsplitting(shiftVal,dShiftVal, switchVal,dSwitchVal);
	shiftVal = switchVal;
	dShiftVal = dSwitchVal;
	break;
      }

//  --------------------------------------------------------------------------
//  END SHIFTING / SPLITTING FUNCTION DEFINITIONS
//  --------------------------------------------------------------------------


      kqq = kq_i * a_j.charge;
      f = kqq*r_1;
      
EXCL
(
      electEnergy -= f * shiftVal;
      if ( m14 ) electEnergy += f * shiftVal * scale14;
)
NOEXCL
(
      electEnergy += f * shiftVal;
)

FULL
(
  EXCL
  (
      fullElectEnergy += f * ( shiftVal - 1.0 );
      if ( m14 ) fullElectEnergy -= f * ( shiftVal - 1.0 ) * scale14;
      BigReal f2 = f * r_1*r_1;
  )
  NOEXCL
  (
      fullElectEnergy -= f * shiftVal;
  )
)
      f *= r_1*(r_1*shiftVal - dShiftVal );

      NOEXCL( const ) BigReal f_elec = f;

EXCL
(
      force_r -= f_elec;
      FULL( fullforce_r += ( f - f2 ); )
      if ( m14 )
      {
	f_elec *= scale14;
	force_r += f_elec;
	FULL( fullforce_r -= ( ( f - f2 ) * scale14 ); )
      }
)
NOEXCL
(
      force_r += f_elec;
      FULL( fullforce_r -= f_elec; )
)

      BigReal r_6 = r_1*r_1*r_1; r_6 *= r_6;
      const BigReal r_12 = r_6*r_6;

      const BigReal A = scaling * lj_pars->A;
      const BigReal B = scaling * lj_pars->B;

      const BigReal AmBterm = (A*r_6 - B) * r_6;

EXCL
(
      vdwEnergy -= switchVal * AmBterm;

      if ( m14 )
      {
	lj_pars = ljTable->table_val(a_i.type, a_j.type, 1);
	const BigReal A = scaling * lj_pars->A;
	const BigReal B = scaling * lj_pars->B;
	const BigReal AmBterm = (A*r_6 - B) * r_6;
	vdwEnergy += switchVal * AmBterm;
	force_r += ( ( switchVal * 6.0 * (A*r_12 + AmBterm) *
			r_1 - AmBterm*dSwitchVal )*r_1 );
      }
)
NOEXCL
(
      vdwEnergy += switchVal * AmBterm;
)



EXCL
(
      force_r -= ( switchVal * 6.0 * (A*r_12 + AmBterm) *
			r_1 - AmBterm*dSwitchVal )*r_1;
)
NOEXCL
(
      force_r += ( switchVal * 6.0 * (A*r_12 + AmBterm) *
			r_1 - AmBterm*dSwitchVal )*r_1;
)

      FAST
      (
      {
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
      }
      )

      FULL
      (
      {
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

NOEXCL
(
    } // for pairlist
  } // for i
  if (pairlist != pairlist_std) delete [] pairlist;
)

  FAST
  (
  reduction[vdwEnergyIndex] += vdwEnergy;
  reduction[electEnergyIndex] += electEnergy;
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

