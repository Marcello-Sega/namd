//-*-c++-*-
/***************************************************************************/
/*                                                                         */
/*              (C) Copyright 1996 The Board of Trustees of the            */
/*                          University of Illinois                         */
/*                           All Rights Reserved                           */
/*									   */
/***************************************************************************/

/***************************************************************************
 * DESCRIPTION: ComputeNonbondedBase.h
 *
 ***************************************************************************/

// Several special cases are defined:
//   NBTYPE: exclusion method (NBPAIR, NBSELF, NBEXCL -- mutually exclusive)
//   FULLELECT full electrostatics calculation?

#ifdef DEFINITION // (
  #include "LJTable.h"
  #include "Molecule.h"
  #include "ComputeNonbondedUtil.h"
#endif // )

// only define this when using hydrogen grouping code.
// don't define this if you want the original code.
#if 1
  #define HGROUPING(X) X
  #define NOHGROUPING(X)
#else
  #define HGROUPING(X)
  #define NOHGROUPING(X) X
#endif

// determining class name
#undef NAME
#undef CLASS
#undef CLASSNAME
#define NAME CLASSNAME(calc)

#undef PAIR
#if NBTYPE == NBPAIR
  #define PAIR(X) X
  #define CLASS ComputeNonbondedPair
  #define CLASSNAME(X) FULLELECTNAME( X ## _pair )
#else
  #define PAIR(X)
#endif

#undef SELF
#if NBTYPE == NBSELF
  #define SELF(X) X
  #define CLASS ComputeNonbondedSelf
  #define CLASSNAME(X) FULLELECTNAME( X ## _self )
#else
  #define SELF(X)
#endif

#undef EXCL
#undef NOEXCL
#if NBTYPE == NBEXCL
  #define EXCL(X) X
  #define CLASS ComputeNonbondedExcl
  #define CLASSNAME(X) FULLELECTNAME( X ## _excl )
  #define NOEXCL(X)
#else
  #define EXCL(X)
  #define NOEXCL(X) X
#endif

#undef FULLELECTNAME
#undef FULL
#undef NOFULL
#ifdef FULLELECT
  #define FULLELECTNAME(X) SPLITTINGNAME( X ## _fullelect )
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
    Force & f_i = params->ff[0][0];
    Force & f_j = params->ff[1][0];
    const AtomProperties & a_i = params->a[0][0];
    const AtomProperties & a_j = params->a[1][0];
    int & m14 = params->m14;
    // used by full electrostatics
    Force & fullf_i = params->fullf[0][0];
    Force & fullf_j = params->fullf[1][0];
  )

  // local variables
  BigReal vdwEnergy = 0;
  BigReal electEnergy = 0;
  BigReal virial = 0;
  FULL
  (
  BigReal fullElectEnergy = 0;
  BigReal fullElectVirial = 0;	// value == fullElectEnergy until the end
  )
  NOEXCL
  (
  // Bringing stuff into local namespace for speed.
  // Probably makes things slower in exclusion mode, though.

  register const BigReal cutoff2 = ComputeNonbondedUtil:: cutoff2;
  HGROUPING
  (
  register const BigReal groupcutoff2 = ComputeNonbondedUtil:: groupcutoff2;
  )
  const BigReal dielectric_1 = ComputeNonbondedUtil:: dielectric_1;

  const LJTable* const ljTable = ComputeNonbondedUtil:: ljTable;
  const Molecule* const mol = ComputeNonbondedUtil:: mol;
  const BigReal scale14 = ComputeNonbondedUtil:: scale14;
  const Real switchOn = ComputeNonbondedUtil:: switchOn;
  const BigReal switchOn2 = ComputeNonbondedUtil:: switchOn2;
  // const BigReal c0 = ComputeNonbondedUtil:: c0;
  const BigReal c1 = ComputeNonbondedUtil:: c1;
  const BigReal c3 = ComputeNonbondedUtil:: c3;
  const BigReal c5 = ComputeNonbondedUtil:: c5;
  const BigReal c6 = ComputeNonbondedUtil:: c6;
  const BigReal d0 = ComputeNonbondedUtil:: d0;
  )

  BigReal kqq;	// initialized later
  BigReal f;	// initialized later
  register BigReal r2;

  // for speedup
  register BigReal tmp_x;
  register BigReal tmp_y;
  register BigReal tmp_z;

EXCL
(
  const BigReal kq_i = COLOUMB * a_i.charge * dielectric_1;
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

  HGROUPING
  (
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
  )

  // for speeding up the for-loop
  const AtomProperties *a_0 = params->a[0];
  const AtomProperties *a_1 = params->a[1];
  const Position *p_0 = params->p[0];
  const Position *p_1 = params->p[1];
  Force *f_0 = params->ff[0];
  Force *f_1 = params->ff[1];
  FULL
    (
    Force *fullf_0 = params->fullf[0];
    Force *fullf_1 = params->fullf[1];
    )

  for ( i = 0; i < i_upper; ++i )
  {
    const AtomProperties &a_i = a_0[i];
    const Position &p_i = p_0[i];
    register const BigReal p_i_x = p_i.x;
    register const BigReal p_i_y = p_i.y;
    register const BigReal p_i_z = p_i.z;

    Force & f_i = f_0[i];
    FULL( Force & fullf_i = fullf_0[i]; )

  HGROUPING
  (
  if (a_i.hydrogenGroupSize) // if hydrogen group parent
    {
    pairlistindex = 0;	// initialize with 0 elements
    pairlistoffset=0;
    int groupfixed = ( a_i.flags & GROUP_FIXED );

    // If patch divisions are not made by hydrogen groups, then
    // hydrogenGroupSize is set to 1 for all atoms.  Thus we can
    // carry on as if we did have groups - only less efficiently.
    // An optimization in this case is to not rebuild the temporary
    // pairlist but to include every atom in it.  This should be a
    // a very minor expense.

    register const Position *p_j = p_1;
    SELF( p_j += i+1; )

    PAIR( j = 0; )
    SELF
      (
      // add all child hydrogens of i
      for( j=i+1; (j<j_upper) && (a_1[j].hydrogenGroupSize == 0); j++)
	{
	pairlist[pairlistindex++] = j;
	p_j++;
	}
      )

    // add remaining atoms to pairlist via hydrogen groups
    register const AtomProperties *pa_j = a_1 + j;
    register BigReal p_j_x = p_j->x;
    register BigReal p_j_y = p_j->y;
    register BigReal p_j_z = p_j->z;
    register int *pli = pairlist + pairlistindex;

    while ( j < j_upper )
	{
	register int hgs = pa_j->hydrogenGroupSize;
	p_j += ( ( j + hgs < j_upper ) ? hgs : 0 );
	r2 = p_i_x - p_j_x;
	r2 *= r2;
	p_j_x = p_j->x;					// preload
	register BigReal t2 = p_i_y - p_j_y;
	r2 += t2 * t2;
	p_j_y = p_j->y;					// preload
	t2 = p_i_z - p_j_z;
	r2 += t2 * t2;
	p_j_z = p_j->z;					// preload
	// use a slightly large cutoff to include hydrogens
	if ( r2 <= groupcutoff2 &&
		! ( groupfixed && (pa_j->flags & GROUP_FIXED) ) )
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
  )

  const BigReal kq_i_u = COLOUMB * a_i.charge * dielectric_1;

    const BigReal kq_i_s = kq_i_u * scale14;
    register const Position *p_j = p_1;

    NOHGROUPING
    (
    SELF
      (
        if ( i + 1 < j_upper ) p_j += i+1;
      )
    )
    HGROUPING
    (
      if ( pairlistoffset < pairlistindex ) p_j += pairlist[pairlistoffset];
    )

    register BigReal p_j_x = p_j->x;
    register BigReal p_j_y = p_j->y;
    register BigReal p_j_z = p_j->z;

    HGROUPING
    (
    for (int k=pairlistoffset; k<pairlistindex; k++)
    {
      j = pairlist[k];
      // don't worry about [k+1] going beyond array since array is 1 too large
      p_j += pairlist[k+1]-j; // preload
    )
    NOHGROUPING
    (

    SELF( j = i+1; )
    PAIR( j = 0; )
    for( ; j<j_upper; j++)
    {
      p_j += ( j + 1 < j_upper );
    )
      register const BigReal p_ij_x = p_i_x - p_j_x;
      p_j_x = p_j->x;					// preload
      register const BigReal p_ij_y = p_i_y - p_j_y;
      p_j_y = p_j->y;					// preload
      register const BigReal p_ij_z = p_i_z - p_j_z;
      p_j_z = p_j->z;					// preload

)

      // common code
      r2 = square(p_ij_x,p_ij_y,p_ij_z);

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
	    reduction[fullElectVirialIndex] += fullElectEnergy;
	  )
	return; )
      }

NOEXCL
(
      BigReal kq_i = kq_i_u;
      const AtomProperties & a_j = a_1[j];

      FULL
      (
        Force & fullf_j = fullf_1[j];
	const BigReal r = sqrt(r2);
        const BigReal r_1 = 1.0/r;
        kqq = kq_i * a_j.charge;
        f = kqq*r_1;
      )
)

      register BigReal force_r = 0.0;		//  force / r
      FULL
      (
      register BigReal fullforce_r = 0.0;	//  fullforce / r
      )

      const LJTable::TableEntry * lj_pars = 
		ljTable->table_val(a_i.type, a_j.type);

      if ( r2 <= lj_pars->exclcut2 )
      {
	NOEXCL
	(
	if ( mol->checkexcl(a_i.id,a_j.id) )  // Inline this by hand.
	{
	  FULL
	  (
	    // Do a quick fix and get out!
	    fullElectEnergy -= f;
	    fullforce_r = -f * r_1*r_1;
	    tmp_x = fullforce_r * p_ij_x;
	    fullf_i.x += tmp_x;
	    fullf_j.x -= tmp_x;
	    tmp_y = fullforce_r * p_ij_y;
	    fullf_i.y += tmp_y;
	    fullf_j.y -= tmp_y;
	    tmp_z = fullforce_r * p_ij_z;
	    fullf_i.z += tmp_z;
	    fullf_j.z -= tmp_z;
	  )
	  continue;  // Must have stored force by now.
	}
	else if ( mol->check14excl(a_i.id,a_j.id) )  // Inline this by hand.
	{
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
	)

        EXCL
        (
	  return;
        )

      }

      NOEXCL
      (
      Force & f_j = f_1[j];
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
	shifting(shiftVal,dShiftVal,r,r2,c5,c6);
	break;
      case SPLIT_C1:
	c1splitting(shiftVal,dShiftVal,r,d0,switchOn);
	break;
      case SPLIT_XPLOR:
	xplorsplitting(shiftVal,dShiftVal, switchVal,dSwitchVal);
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

      const BigReal &A = lj_pars->A;
      const BigReal &B = lj_pars->B;

      const BigReal AmBterm = (A*r_6 - B) * r_6;

EXCL
(
      vdwEnergy -= switchVal * AmBterm;

      if ( m14 )
      {
	lj_pars = ljTable->table_val(a_i.type, a_j.type, 1);
	const BigReal &A = lj_pars->A;
	const BigReal &B = lj_pars->B;
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

      virial += force_r * r2;

      tmp_x = force_r * p_ij_x;
      f_i.x += tmp_x;
      f_j.x -= tmp_x;
      tmp_y = force_r * p_ij_y;
      f_i.y += tmp_y;
      f_j.y -= tmp_y;
      tmp_z = force_r * p_ij_z;
      f_i.z += tmp_z;
      f_j.z -= tmp_z;

      FULL
      (
      fullElectVirial = fullElectEnergy + fullforce_r * r2;

      tmp_x = fullforce_r * p_ij_x;
      fullf_i.x += tmp_x;
      fullf_j.x -= tmp_x;
      tmp_y = fullforce_r * p_ij_y;
      fullf_i.y += tmp_y;
      fullf_j.y -= tmp_y;
      tmp_z = fullforce_r * p_ij_z;
      fullf_i.z += tmp_z;
      fullf_j.z -= tmp_z;
      )

NOEXCL
(
    } // for pairlist
  } // for i
  HGROUPING( if (pairlist != pairlist_std) delete [] pairlist; )
)

  reduction[vdwEnergyIndex] += vdwEnergy;
  reduction[electEnergyIndex] += electEnergy;
  reduction[virialIndex] += virial;
  FULL
  (
  reduction[fullElectEnergyIndex] += fullElectEnergy;
  reduction[fullElectVirialIndex] += fullElectVirial;
  )
}

/***************************************************************************
 * RCS INFORMATION:
 *
 *	$RCSfile: ComputeNonbondedBase.h,v $
 *	$Author: jim $	$Locker:  $		$State: Exp $
 *	$Revision: 1.1028 $	$Date: 1997/09/19 08:55:30 $
 *
 ***************************************************************************
 * REVISION HISTORY:
 *
 * $Log: ComputeNonbondedBase.h,v $
 * Revision 1.1028  1997/09/19 08:55:30  jim
 * Added rudimentary but relatively efficient fixed atoms.  New options
 * are fixedatoms, fixedatomsfile, and fixedatomscol (nonzero means fixed).
 * Energies will be affected, although this can be fixed with a little work.
 *
 * Revision 1.1027  1997/09/19 05:17:42  jim
 * Cleaned up and tweaked hydrogen-group based temporary pairlist
 * generation for roughly a 6% performance improvement.
 *
 * Revision 1.1026  1997/09/19 02:59:56  jim
 * Added caching of groupcutoff2.
 *
 * Revision 1.1025  1997/09/05 20:14:27  jim
 * Small fixes.
 *
 * Revision 1.1024  1997/07/30 21:23:22  jim
 * More possible bugs (memory references).
 *
 * Revision 1.1023  1997/07/30 20:51:47  jim
 * Probable bug fix - added sentinal to end of parilist.
 *
 * Revision 1.1022  1997/06/05 20:19:41  nealk
 * Minor modifications for readability and very minor speedup.
 *
 * Revision 1.1021  1997/06/04 20:13:50  nealk
 * Modified to simplify macros.
 *
 * Revision 1.1020  1997/05/29 19:14:00  nealk
 * Removed some array indexing for minor speed improvement.
 *
 * Revision 1.1019  1997/05/23 19:29:38  nealk
 * Removed more macros.
 *
 * Revision 1.1018  1997/05/20 15:49:08  nealk
 * Pair, Self, and Excl not use the same parameters!
 *
 * Revision 1.1016  1997/05/13 18:30:45  nealk
 * Removed ComputeNonbondedHack.h!
 * Reduced a lot of code in Util and Base.
 * ComputeNonbondedBase.h now only contains the function definitions.
 * The only heavy macro areas are in Util.C (determining which Base.h to define)
 * and Base.h (where the functions are defined).
 *
 * Revision 1.1015  1997/05/12 18:45:21  nealk
 * Minor coding changes (looks nicer).
 *
 * Revision 1.1014  1997/05/09 18:24:22  nealk
 * 1. Added hydrogen grouping code to improve performance in ComputeNonbondedBase
 *    CODE ONLY WORKS WITH HYDROGEN GROUPING!
 * 2. Increased the hydrogen group cutoff side from 2A to 2.5A -- 2A gave
 *    fractionally different values after 100 iterations.  2.5A gives same numbers.
 * 3. Made migration by hydrogen grouping the default in SimParameters.
 *
 * Revision 1.1013  1997/05/05 16:38:57  nealk
 * Corrected cutoff value used with hydrogen grouping.  (groupcutoff2)
 *
 * Revision 1.1012  1997/05/05 15:28:24  nealk
 * Added water-water specific code to NonbondedBase.  The cutoff for the temp
 * pairlist is currently disabled.
 *
 * Revision 1.1011  1997/03/17 03:55:23  jim
 * Reordered final force store operations for better work/memory interleaving.
 *
 * Revision 1.1010  1997/03/17 03:44:14  jim
 * Rearranged final force store for better memory access (I hope).
 *
 * Revision 1.1009  1997/03/17 03:15:01  jim
 * Added virial calculation.
 *
 * Revision 1.1008  1997/03/17 02:52:13  jim
 * Did cleanups and speedups in preparation for virial calculation.
 * Fixed minor bug which resulted in incorrect energies for excluded
 * pairs which were outside of the cutoff radius (VERY rare).
 *
 * Revision 1.1007  1997/03/14 23:18:08  jim
 * Implemented C1 splitting for long-range electrostatics.
 * Energies on sub-timesteps are incorrect.  (Aren't they always?)
 *
 * Revision 1.1006  1997/03/14 06:44:53  jim
 * First working versions of full electrostatics splitting functions.
 *
 * Revision 1.1005  1997/03/10 00:49:46  jim
 * Eliminated constant copying in exclusion mode, hopefully saves time.
 *
 * Revision 1.1004  1997/02/28 04:47:02  jim
 * Full electrostatics now works with fulldirect on one node.
 *
 * Revision 1.1003  1997/02/21 20:45:11  jim
 * Eliminated multiple function for switching and modified 1-4 interactions.
 * Now assumes a switching function, but parameters are such that nothing
 * happens, same for modified 1-4.  Slight penalty for rare simulations
 * in which these features are not used, but otherwise no loss and
 * simplifies code.
 *
 * Revision 1.1002  1997/02/13 23:17:15  ari
 * Fixed a final bug in AtomMigration - numatoms in ComputePatchPair.C not
 * set correctly in atomUpdate()
 *
 * Revision 1.1001  1997/02/10 08:30:28  jim
 * Now handles periodic boundaries correctly (I forgot this when I was
 * changing Bonds, Angles, etc.) but a bit of a hack, needs cleaning up.
 *
 * Revision 1.1000  1997/02/06 15:58:06  ari
 * Resetting CVS to merge branches back into the main trunk.
 * We will stick to main trunk development as suggested by CVS manual.
 * We will set up tags to track fixed points of development/release
 * as suggested by CVS manual - all praise the CVS manual.
 *
 * Revision 1.778  1997/01/28 00:30:18  ari
 * internal release uplevel to 1.778
 *
 * Revision 1.777.2.5  1997/01/27 22:45:04  ari
 * Basic Atom Migration Code added.
 * Added correct magic first line to .h files for xemacs to go to C++ mode.
 * Compiles and runs without migration turned on.
 *
 * Revision 1.777.2.4  1997/01/23 05:13:14  jim
 * Added pre-fetching of positions for cutoff check loop.
 *
 * Revision 1.777.2.3  1997/01/23 04:20:50  jim
 * Converted cutoff check to use register variables.
 *
 * Revision 1.777.2.2  1997/01/22 21:42:11  jim
 * Larger patches, no two-away computes, small tweak to inner loop.
 *
 * Revision 1.777.2.1  1997/01/17 20:48:01  jim
 * Fixed Lennard-Jones error, now runs with identical energies for e timesteps.
 *
 * Revision 1.777  1997/01/17 19:35:52  ari
 * Internal CVS leveling release.  Start development code work
 * at 1.777.1.1.
 *
 * Revision 1.11  1997/01/17 19:27:00  jim
 * Now energies are reported correctly, electrostatic energy agrees with
 * namd 1.X but Lennard-Jones is still off.
 *
 * Revision 1.10  1997/01/16 19:59:56  jim
 * Added reduction calls to ComputeNonbondedSelf and ...Pair.
 * Also moved some code from ...Excl to ...Util.
 *
 * Revision 1.9  1996/12/04 17:16:32  jim
 * ComputeNonbondedUtil::select() now caches simulation parameters
 *
 * Revision 1.8  1996/12/03 21:05:09  jim
 * added support for exclusion correction computes
 *
 * Revision 1.7  1996/11/30 20:30:36  jim
 * turned off some debugging, switched to DebugM()
 *
 * Revision 1.6  1996/11/23 23:00:30  jim
 * added debug message - energy output
 *
 * Revision 1.5  1996/11/20 23:16:39  jim
 * first compiling version of generic nonbonded function
 *
 * Revision 1.4  1996/11/08 02:37:12  jim
 * split into two files to hide some preprocessor code from user
 *
 *
 ***************************************************************************/

