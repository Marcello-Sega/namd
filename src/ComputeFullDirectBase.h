/***************************************************************************/
/*       (C) Copyright 1996,1997 The Board of Trustees of the              */
/*                          University of Illinois                         */
/*                           All Rights Reserved                           */
/***************************************************************************/
/***************************************************************************
 * DESCRIPTION:  Function body for calc_fulldirect.
 *
 ***************************************************************************/

{
  const BigReal coloumb = COLOUMB * ComputeNonbondedUtil::dielectric_1;
  BigReal *dp1 = data1;
  BigReal *rp1 = results1;
  int j_begin = 0;
  register BigReal electEnergy = 0.;

#ifdef FULLDIRECT_PERIODIC
  BigReal a1 = lattice->a();  BigReal b1 = ( a1 ? 1. / a1 : 0 );
  BigReal a2 = lattice->b();  BigReal b2 = ( a2 ? 1. / a2 : 0 );
  BigReal a3 = lattice->c();  BigReal b3 = ( a3 ? 1. / a3 : 0 );
#endif

  for(int i=0; i<n1; ++i)
  {
    register BigReal p_i_x = *(dp1++);
    register BigReal p_i_y = *(dp1++);
    register BigReal p_i_z = *(dp1++);
    register BigReal kq_i = coloumb * *(dp1++);
    register BigReal f_i_x = 0.;
    register BigReal f_i_y = 0.;
    register BigReal f_i_z = 0.;
    if ( selfmode )
    {
      ++j_begin; data2 += 4; results2 += 3;
    }
    register BigReal *dp2 = data2;
    register BigReal *rp2 = results2;
    register int n2c = n2;
    register int j;
    for( j = j_begin; j<n2c; ++j)
    {
      register BigReal p_ij_x = p_i_x - *(dp2++);
      register BigReal p_ij_y = p_i_y - *(dp2++);
      register BigReal p_ij_z = p_i_z - *(dp2++);

#ifdef FULLDIRECT_PERIODIC
      p_ij_x -= a1 * rint( b1 * p_ij_x );
      p_ij_y -= a2 * rint( b2 * p_ij_y );
      p_ij_z -= a3 * rint( b3 * p_ij_z );
#endif

      register BigReal r_1;
      r_1 = 1./sqrt(p_ij_x * p_ij_x + p_ij_y * p_ij_y + p_ij_z * p_ij_z);
      register BigReal f = *(dp2++) * kq_i * r_1;
      electEnergy += f;
      f *= r_1*r_1;
      p_ij_x *= f;
      p_ij_y *= f;
      p_ij_z *= f;
      f_i_x += p_ij_x;
      f_i_y += p_ij_y;
      f_i_z += p_ij_z;
      *(rp2++) -= p_ij_x;
      *(rp2++) -= p_ij_y;
      *(rp2++) -= p_ij_z;
    }
    *(rp1++) += f_i_x;
    *(rp1++) += f_i_y;
    *(rp1++) += f_i_z;
  }

  return electEnergy;
}

/***************************************************************************
 * RCS INFORMATION:
 *
 *	$RCSfile $
 *	$Author $	$Locker:  $		$State: Exp $
 *	$Revision: 1.1 $	$Date: 1998/03/30 21:01:16 $
 *
 ***************************************************************************
 * REVISION HISTORY:
 *
 * $Log: ComputeFullDirectBase.h,v $
 * Revision 1.1  1998/03/30 21:01:16  jim
 * Added nearest-image support for periodic boundary conditions to full direct.
 *
 *
 ***************************************************************************/
