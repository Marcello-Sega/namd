/***************************************************************************/
/*                                                                         */
/*              (C) Copyright 1996 The Board of Trustees of the            */
/*                          University of Illinois                         */
/*                           All Rights Reserved                           */
/*									   */
/***************************************************************************/

/***************************************************************************
 * DESCRIPTION: ComputeNonbondedPair.C
 *
 ***************************************************************************/

#include "ComputeNonbondedPair.h"
#include "LJTable.h"
#include "Node.h"
#include "SimParameters.h"
#include "Molecule.h"

void ComputeNonbondedPair::doForce(Position* p[2],
                               Force* f[2],
                               AtomProperties* a[2])
{
    CPrintf("ComputeNonbondedPair::doForce() - Eval was sent\n");
    CPrintf(" %d patch 1 atoms and %d patch 2 atoms\n", numAtoms[0], numAtoms[1] );

   Vector f_elec;                       //  Electrostatic force vector
   Vector f_vdw;                        //  Vdw force vector
					//  same as vector from atom2 to atom1
   BigReal r;                           //  Distance between current pair of
                                        //  atoms
   BigReal r2;                          //  r squared
   BigReal r_1;                         //  1/r
   BigReal r_6, r_12;                   //  1/r^6 and 1/r^12
   BigReal switchVal;                   //  Value of switching function
   BigReal shiftVal;                    //  Value of shifting function
   BigReal dShiftVal;                   //  d(shiftVal)/dr
   BigReal dSwitchVal;                  //  d(switchVal)/dr

  BigReal electEnergy = 0;
  BigReal vdwEnergy = 0;

  const int excludeFlag  = Node::Object()->simParameters->exclude;
  const BigReal scale14 = ( excludeFlag == SCALED14 ?
                   Node::Object()->simParameters->scale14 : 1 );

  const Real cutoff = Node::Object()->simParameters->cutoff;
  const BigReal cutoff2 = cutoff*cutoff;
  const Real switchOn = Node::Object()->simParameters->switchingDist;
  const BigReal switchOn2 = switchOn*switchOn;

  const LJTable* const ljTable = LJTable::Instance();
  const Molecule* const mol = Node::Object()->molecule;
  const BigReal dielectric_1 = 1/Node::Object()->simParameters->dielectric;

  const BigReal c0 = 1/(cutoff2-switchOn2);
  const BigReal c1 = c0*c0*c0;
  const BigReal c3 = c1 * 4;
  const BigReal c5 = 1/cutoff2;
  const BigReal c6 = -4 * c5;

  for ( int i = 0; i < numAtoms[0]; ++i )
  {
    const Position & p_i = p[0][i];
    const AtomProperties & a_i = a[0][i];

    Force & f_i = f[0][i];

    const BigReal kq_i = COLOUMB * a_i.charge * dielectric_1;

    for ( int j = 0; j < numAtoms[1]; ++j )
    {
      const Position & p_j = p[1][j];

      f_vdw = p_i - p_j;

      if ( ( r2 = f_vdw.length2() ) > cutoff2 ) continue;

      const AtomProperties & a_j = a[1][j];

      const LJTable::TableEntry * lj_pars = 
		ljTable->table_val(a_i.type, a_j.type);

      BigReal scaleFactor = 1.;

      if ( r2 <= lj_pars->exclcut2 )
      {
	if ( mol->checkexcl(a_i.id,a_j.id) ) continue;
	else if ( ( excludeFlag == SCALED14 ) &&
				mol->check14excl(a_i.id,a_j.id) )
	{
	  lj_pars = ljTable->table_val(a_i.type, a_j.type, 1);
	  scaleFactor = scale14;
	}
      }

      Force & f_j = f[1][j];

      r = sqrt(r2);
      r_1 = 1/r;

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
 *	$RCSfile: ComputeNonbondedPair.C,v $
 *	$Author: jim $	$Locker:  $		$State: Exp $
 *	$Revision: 1.6 $	$Date: 1996/11/05 21:12:12 $
 *
 ***************************************************************************
 * REVISION HISTORY:
 *
 * $Log: ComputeNonbondedPair.C,v $
 * Revision 1.6  1996/11/05 21:12:12  jim
 * fixed modified pairs
 *
 * Revision 1.5  1996/11/05 05:08:56  jim
 * added nonbonded compute code for one case ( no ifs )
 *
 * Revision 1.4  1996/10/31 21:43:29  jim
 * First incarnation as ...Pair
 *
 * Revision 1.3  1996/10/30 01:16:32  jim
 * added AtomProperties structure in Patch plus boxes, passing, etc.
 *
 * Revision 1.2  1996/10/30 00:16:16  jim
 * Removed PositionArray usage.
 *
 * Revision 1.1  1996/10/29 23:55:54  jim
 * Initial revision
 *
 *
 ***************************************************************************/

