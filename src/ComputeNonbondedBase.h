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
//   EXCLUDE14 modified 1-4 parameters?
//   SWITCHING switching function?
//   FULLELECT full electrostatics?

// Check mutual exclusivity
#ifdef NBPAIR
#ifdef NBSELF
#error Error!
#endif
#ifdef NBEXCL
#error Error!
#endif
#endif
#ifdef NBSELF
#ifdef NBEXCL
#error Error!
#endif
#endif

#if defined NBPAIR
#define PLEN [2]
#define I_SUB 0][i
#define I_LOWER 0
#define I_UPPER numAtoms[0]
#define J_SUB 1][j
#define J_LOWER 0
#define J_UPPER numAtoms[1]
#elif defined NBSELF
#define PLEN
#define I_SUB [i]
#define I_LOWER 0
#define I_UPPER numAtoms
#define J_SUB [j]
#define J_LOWER i + 1
#define J_UPPER numAtoms
#endif

(Position* p PLEN, Force* f PLEN, AtomProperties* a PLEN)
{
  Vector f_elec;			//  Electrostatic force vector
  Vector f_vdw;				//  Vdw force vector
					//  same as vector from atom2 to atom1
  BigReal r;				//  Distance between current pair of
					//  atoms
  BigReal r2;				//  r squared
  BigReal r_1;				//  1/r
  BigReal r_6, r_12;			//  1/r^6 and 1/r^12
   
  BigReal electEnergy = 0;
  BigReal vdwEnergy = 0;

  const Real cutoff = Node::Object()->simParameters->cutoff;
  const BigReal cutoff2 = cutoff*cutoff;

  const LJTable* const ljTable = LJTable::Instance();
  const Molecule* const mol = Node::Object()->molecule;
  const BigReal dielectric_1 = 1/Node::Object()->simParameters->dielectric;

#ifdef EXCLUDE14
  const BigReal scale14 = Node::Object()->simParameters->scale14;
#endif

#ifdef SWITCHING
  const Real switchOn = Node::Object()->simParameters->switchingDist;
  const BigReal switchOn2 = switchOn*switchOn;
  const BigReal c0 = 1/(cutoff2-switchOn2);
  const BigReal c1 = c0*c0*c0;
  const BigReal c3 = c1 * 4;
  const BigReal c5 = 1/cutoff2;
  const BigReal c6 = -4 * c5;
  BigReal switchVal;                   //  Value of switching function
  BigReal shiftVal;                    //  Value of shifting function
  BigReal dShiftVal;                   //  d(shiftVal)/dr
  BigReal dSwitchVal;                  //  d(switchVal)/dr
#else
#define switchVal 1
#define shiftVal 1
#define dShiftVal 0
#define dSwitchVal 0
#endif

  for ( int i = I_LOWER; i < I_UPPER; ++i )
  {
    const Position & p_i = p[I_SUB];
    const AtomProperties & a_i = a[I_SUB];

    Force & f_i = f[I_SUB];

    const BigReal kq_i = COLOUMB * a_i.charge * dielectric_1;

    for ( int j = J_LOWER; j < J_UPPER; ++j )
    {
      const Position & p_j = p[J_SUB];

      f_vdw = p_i - p_j;

      if ( ( r2 = f_vdw.length2() ) > cutoff2 ) continue;

      const AtomProperties & a_j = a[J_SUB];

      const LJTable::TableEntry * lj_pars = 
		ljTable->table_val(a_i.type, a_j.type);
      BigReal scaleFactor = 1.;

      if ( r2 <= lj_pars->exclcut2 )
      {
	if ( mol->checkexcl(a_i.id,a_j.id) ) continue;
	
#ifdef EXCLUDE14
	else if ( mol->check14excl(a_i.id,a_j.id) )
	{
	  lj_pars = ljTable->table_val(a_i.type, a_j.type, 1);
	  scaleFactor = scale14;
	}
#endif

      }

      Force & f_j = f[J_SUB];

      r = sqrt(r2);
      r_1 = 1/r;

#ifdef SWITCHING
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
#endif

      const BigReal kqq = kq_i * a_j.charge * scaleFactor;

      BigReal f = kqq*r_1;
      
      electEnergy += f*shiftVal;
      
      f *= r_1*(shiftVal*r_1 - dShiftVal);

      f_elec = f_vdw*f;

      f_i += f_elec;
      f_j -= f_elec;

      r_6 = r_1*r_1*r_1;
      r_6 *= r_6;
      r_12 = r_6*r_6;

      const BigReal A = lj_pars->A;
      const BigReal B = lj_pars->B;

      const BigReal AmBterm = (A*r_6 - B) * r_6;

      vdwEnergy += switchVal*AmBterm;
      
      f_vdw *= (switchVal* 6 * (A*r_12 + AmBterm)*r_1-AmBterm*dSwitchVal)*r_1;

      f_i += f_vdw;
      f_j -= f_vdw;
    }
  }
}

/***************************************************************************
 * RCS INFORMATION:
 *
 *	$RCSfile: ComputeNonbondedBase.h,v $
 *	$Author: jim $	$Locker:  $		$State: Exp $
 *	$Revision: 1.1 $	$Date: 1996/11/06 21:17:11 $
 *
 ***************************************************************************
 * REVISION HISTORY:
 *
 ***************************************************************************/

