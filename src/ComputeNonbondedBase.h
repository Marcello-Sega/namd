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
//   NBPAIR, NBSELF, NBEXCL switch environment (mutually exclusive)
//   MODIFY14 modified 1-4 parameters?
//   SWITCHING switching function?

#include "ComputeNonbondedHack.h"

#ifdef DEFINITION

#include "LJTable.h"
#include "Molecule.h"

#define MIN_DEBUG_LEVEL 5
#define DEBUGM
#include "Debug.h"

#endif

#if defined DECLARATION || defined DEFINITION
DECL( static ) void NODECL( ComputeNonbondedUtil :: ) NAME
NOEXCL
(
(Position* p PLEN, Force* f PLEN, AtomProperties* a PLEN,
 int numAtoms PLEN, BigReal *reduction)
)
EXCL
(
(const Position & p_i, const Position & p_j,
 Force & f_i, Force & f_j,
 const AtomProperties & a_i, const AtomProperties & a_j,
 int m14, BigReal *reduction)
)
DECL( ; )
#endif
#ifdef DEFINITION
{
  BigReal electEnergy = 0;
  BigReal vdwEnergy = 0;

  const BigReal cutoff2 = ComputeNonbondedUtil:: cutoff2;
  const BigReal dielectric_1 = ComputeNonbondedUtil:: dielectric_1;

  const LJTable* const ljTable = ComputeNonbondedUtil:: ljTable;
  const Molecule* const mol = ComputeNonbondedUtil:: mol;

M14
(
  const BigReal scale14 = ComputeNonbondedUtil:: scale14;
)

SW
(
  const Real switchOn = ComputeNonbondedUtil:: switchOn;
  const BigReal switchOn2 = ComputeNonbondedUtil:: switchOn2;
  const BigReal c0 = ComputeNonbondedUtil:: c0;
  const BigReal c1 = ComputeNonbondedUtil:: c1;
  const BigReal c3 = ComputeNonbondedUtil:: c3;
  const BigReal c5 = ComputeNonbondedUtil:: c5;
  const BigReal c6 = ComputeNonbondedUtil:: c6;
)

NOEXCL
(
  for ( int i = I_LOWER; i < I_UPPER; ++i )
  {
    const Position & p_i = p[I_SUB];
    const AtomProperties & a_i = a[I_SUB];

    Force & f_i = f[I_SUB];
)

    const BigReal NOEXCL( M14( kq_i_u ) NOM14( kq_i ) ) EXCL( kq_i ) =
    			COLOUMB * a_i.charge * dielectric_1;

NOEXCL
(
    M14( const BigReal kq_i_s = kq_i_u * scale14; )

    for ( int j = J_LOWER; j < J_UPPER; ++j )
    {
      const Position & p_j = p[J_SUB];
)

      Vector f_vdw = p_i - p_j;

      const BigReal r2 = f_vdw.length2();

      if ( r2 > cutoff2 ) NOEXCL( continue ) EXCL( return );

NOEXCL
(
      M14 ( BigReal kq_i = kq_i_u; )

      const AtomProperties & a_j = a[J_SUB];
)

      const LJTable::TableEntry * NOM14( const ) lj_pars = 
		ljTable->table_val(a_i.type, a_j.type);

      if ( r2 <= lj_pars->exclcut2 )
      {
NOEXCL
(
	if ( mol->checkexcl(a_i.id,a_j.id) ) continue;
	
M14
(
	else if ( mol->check14excl(a_i.id,a_j.id) )
	{
	  lj_pars = ljTable->table_val(a_i.type, a_j.type, 1);
	  kq_i = kq_i_s;
	}
)

)
EXCL
(
	return;
)
      }

NOEXCL
(
      Force & f_j = f[J_SUB];
)

      const BigReal r = sqrt(r2);
      const BigReal r_1 = 1/r;

SW
(
      BigReal switchVal;
      BigReal shiftVal;
      BigReal dSwitchVal;
      BigReal dShiftVal;

      if (r > switchOn)
      {
	const BigReal c2 = cutoff2-r2;
	const BigReal c4 = c2*(cutoff2+2*r2-3*switchOn2);
	switchVal = c2*c4*c1;
	dSwitchVal = c3*r*(c2*c2-c4);
      }
      else
      {
	switchVal = 1;
	dSwitchVal = 0;
      }

      shiftVal = 1 - r2*c5;
      dShiftVal = c6*shiftVal*r;
      shiftVal *= shiftVal;
)

      const BigReal kqq = kq_i * a_j.charge;

      BigReal f = kqq*r_1;
      
EXCL
(
      electEnergy -= f SW( * shiftVal );
M14
(
      if ( m14 ) electEnergy += f SW( * shiftVal ) * scale14;
)
)
NOEXCL
(
      electEnergy += f SW( * shiftVal );
)
      
      f *= r_1*( SW( shiftVal * ) r_1 SW( - dShiftVal ) );

      NOEXCL( const ) Vector f_elec = f_vdw*f;

EXCL
(
      f_i -= f_elec;
      f_j += f_elec;
M14
(
      if ( m14 )
      {
	f_elec *= scale14;
	f_i += f_elec;
	f_j -= f_elec;
      }
)
)
NOEXCL
(
      f_i += f_elec;
      f_j -= f_elec;
)

      BigReal r_6 = r_1*r_1*r_1; r_6 *= r_6;
      const BigReal r_12 = r_6*r_6;

      const BigReal A = lj_pars->A;
      const BigReal B = lj_pars->B;

      const BigReal AmBterm = (A*r_6 - B) * r_6;

EXCL
(
      vdwEnergy -= SW( switchVal * ) AmBterm;
M14
(
      if ( m14 )
      {
	lj_pars = ljTable->table_val(a_i.type, a_j.type, 1);
	const BigReal Ax = lj_pars->A;
	const BigReal Bx = lj_pars->B;
	const BigReal AmBtermx = (A*r_6 - B) * r_6;
	vdwEnergy += SW( switchVal * ) AmBtermx;
	const Vector f_vdwx = f_vdw * ( ( SW( switchVal * ) 6 * (A*r_12 + AmBterm) *
			r_1 SW( - AmBterm*dSwitchVal ) )*r_1 );
	f_i += f_vdwx;
	f_j -= f_vdwx;
      }
)
)
NOEXCL
(
      vdwEnergy += SW( switchVal * ) AmBterm;
)


      f_vdw *= ( SW( switchVal * ) 6 * (A*r_12 + AmBterm) *
			r_1 SW( - AmBterm*dSwitchVal ) )*r_1;

EXCL
(
      f_i -= f_vdw;
      f_j += f_vdw;
)
NOEXCL
(
      f_i += f_vdw;
      f_j -= f_vdw;
)

NOEXCL
(
    }
  }
)

  reduction[electEnergyIndex] += electEnergy;
  reduction[vdwEnergyIndex] += vdwEnergy;

  DebugM(3, "Nonbonded computation results: electEnergy = " 
    << electEnergy << " vdwEnergy = " << vdwEnergy << endl);
}
#endif

/***************************************************************************
 * RCS INFORMATION:
 *
 *	$RCSfile: ComputeNonbondedBase.h,v $
 *	$Author: jim $	$Locker:  $		$State: Exp $
 *	$Revision: 1.11 $	$Date: 1997/01/17 19:27:00 $
 *
 ***************************************************************************
 * REVISION HISTORY:
 *
 * $Log: ComputeNonbondedBase.h,v $
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

