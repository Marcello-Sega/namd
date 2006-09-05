/**
***  Copyright (c) 1995, 1996, 1997, 1998, 1999, 2000 by
***  The Board of Trustees of the University of Illinois.
***  All rights reserved.
**/

EXCLUDED( FAST( foo bar ) )
EXCLUDED( MODIFIED( foo bar ) )
EXCLUDED( NORMAL( foo bar ) )
NORMAL( MODIFIED( foo bar ) )


#ifdef ARCH_POWERPC
     __alignx(16, table_four);
     __alignx(16, p_1);
#pragma unroll(1)
#endif

#pragma ivdep
    for (k=0; k<npairi; ++k) {      

      int table_i = (r2iilist[2*k] >> 14) + r2_delta_expc;  // table_i >= 0 
      const int j = pairlisti[k];
      register const CompAtom *p_j = p_1 + j;
      
      BigReal diffa = r2list[k] - r2_table[table_i];
      const BigReal* const table_four_i = table_four + 16*table_i;

      FAST(
      const LJTable::TableEntry * lj_pars = 
              lj_row + 2 * vdwtype_array[j]  MODIFIED(+ 1);
      )
	
#ifdef ARCH_POWERPC
#if ( SHORT( FAST( 1+ ) ) 0 ) 
#pragma disjoint (*table_four_i, *f_1)
#pragma disjoint (*p_j,          *f_1)
#endif
	
#if ( FULL( 1+ ) 0 )
#pragma disjoint (*table_four_i, *fullf_1)
#pragma disjoint (*p_j,          *fullf_1)
#ifdef f_1
#pragma disjoint (*f_1    , *fullf_1)
#pragma disjoint (*fullf_1, *f_1)
#endif
#endif

      __alignx(16, table_four_i);
      __alignx(16, p_j);
      FAST (
      __alignx(16, lj_pars);
      )
#endif

      /*
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
      */

      BigReal kqq = kq_i * p_j->charge;

      FEP(
      int jfep_type = p_j->partition;
      BigReal lambda_pair = lambda_table_i[2*jfep_type];
      BigReal d_lambda_pair = lambda_table_i[2*jfep_type+1];
      )

      LES( BigReal lambda_pair = lambda_table_i[p_j->partition]; )

#if ( FAST(1+) 0 )
      const BigReal A = scaling * lj_pars->A;
      const BigReal B = scaling * lj_pars->B;

      BigReal vdw_d = A * table_four_i[0] - B * table_four_i[2];
      BigReal vdw_c = A * table_four_i[1] - B * table_four_i[3];
      BigReal vdw_b = A * table_four_i[4] - B * table_four_i[6];
      BigReal vdw_a = A * table_four_i[5] - B * table_four_i[7];

      ENERGY(
      register BigReal vdw_val =
        ( ( diffa * vdw_d * (1/6.)+ vdw_c * (1/4.)) * diffa + vdw_b *(1/2.)) * diffa + vdw_a;
      vdwEnergy -= LAM(lambda_pair *) vdw_val;
      FEP( vdwEnergy_s -= d_lambda_pair * vdw_val; )
      )
#endif // FAST

      register const BigReal p_ij_x = p_i_x - p_j->position.x;
      register const BigReal p_ij_y = p_i_y - p_j->position.y;
      register const BigReal p_ij_z = p_i_z - p_j->position.z;

#if ( FAST(1+) 0 )
      INT( 
      register BigReal vdw_dir =
      ( diffa * vdw_d + vdw_c ) * diffa + vdw_b;
      //BigReal force_r =  LAM(lambda_pair *) vdw_dir;
      reduction[pairVDWForceIndex_X] += force_sign * vdw_dir * p_ij_x;
      reduction[pairVDWForceIndex_Y] += force_sign * vdw_dir * p_ij_y;
      reduction[pairVDWForceIndex_Z] += force_sign * vdw_dir * p_ij_z;
      )

#if ( SHORT(1+) 0 )
      NORMAL(
      BigReal fast_d = kqq * table_four_i[8];
      BigReal fast_c = kqq * table_four_i[9];
      BigReal fast_b = kqq * table_four_i[10];
      BigReal fast_a = kqq * table_four_i[11];
      )
      MODIFIED(
      BigReal modfckqq = (1.0-modf_mod) * kqq;
      BigReal fast_d = modfckqq * table_four_i[8];
      BigReal fast_c = modfckqq * table_four_i[9];
      BigReal fast_b = modfckqq * table_four_i[10];
      BigReal fast_a = modfckqq * table_four_i[11];
      )

      {
      ENERGY(
	     register BigReal fast_val =
	( ( diffa * fast_d * (1/6.)+ fast_c * (1/4.)) * diffa + fast_b *(1/2.)) * diffa + fast_a;
      electEnergy -=  LAM(lambda_pair *) fast_val;
      FEP( electEnergy_s -=  d_lambda_pair * fast_val; )
      )

      INT(
      register BigReal fast_dir =
      ( diffa * fast_d + fast_c ) * diffa + fast_b;
      // force_r -= -1.0 * LAM(lambda_pair *) fast_dir;
      reduction[pairElectForceIndex_X] +=  force_sign * fast_dir * p_ij_x;
      reduction[pairElectForceIndex_Y] +=  force_sign * fast_dir * p_ij_y;
      reduction[pairElectForceIndex_Z] +=  force_sign * fast_dir * p_ij_z;
      )
      }

      fast_d += vdw_d;
      fast_c += vdw_c;
      fast_b += vdw_b;
      fast_a += vdw_a;
      register BigReal fast_dir =
	( diffa * fast_d + fast_c ) * diffa + fast_b;
      BigReal force_r =  LAM(lambda_pair *) fast_dir;
      register BigReal tmp_x = force_r * p_ij_x;
      PAIR( virial_xx += tmp_x * p_ij_x; )
      PAIR( virial_xy += tmp_x * p_ij_y; )
      PAIR( virial_xz += tmp_x * p_ij_z; )

      f_i_x += tmp_x;
      f_1[j].x -= tmp_x;

      register BigReal tmp_y = force_r * p_ij_y;
      PAIR( virial_yy += tmp_y * p_ij_y; )
      PAIR( virial_yz += tmp_y * p_ij_z; )
      f_i_y += tmp_y;
      f_1[j].y -= tmp_y;
      
      register BigReal tmp_z = force_r * p_ij_z;
      PAIR( virial_zz += tmp_z * p_ij_z; )
      f_i_z += tmp_z;
      f_1[j].z -= tmp_z;

      PPROF(
        const BigReal p_j_z = p_j->position.z;
        int n2 = (int)floor((p_j_z-pressureProfileMin)*invThickness);
        pp_clamp(n2, pressureProfileSlabs);
        int p_j_partition = p_j->partition;

        pp_reduction(pressureProfileSlabs, n1, n2, 
                     p_i_partition, p_j_partition, pressureProfileAtomTypes,
                     tmp_x*p_ij_x, tmp_y * p_ij_y, tmp_z*p_ij_z,
                     pressureProfileReduction);

      )

#endif // SHORT
#endif // FAST

      FULL(
      BigReal slow_d = table_four_i[8 SHORT(+ 4)];
      BigReal slow_c = table_four_i[9 SHORT(+ 4)];
      BigReal slow_b = table_four_i[10 SHORT(+ 4)];
      BigReal slow_a = table_four_i[11 SHORT(+ 4)];
      EXCLUDED(
      SHORT(
      const BigReal* const slow_i = slow_table + 4*table_i;
      slow_a +=    slow_i[0];
      slow_b += 2.*slow_i[1];
      slow_c += 4.*slow_i[2];
      slow_d += 6.*slow_i[3];
      )
      NOSHORT(
      slow_d -= table_four_i[12];
      slow_c -= table_four_i[13];
      slow_b -= table_four_i[14];
      slow_a -= table_four_i[15];
      )
      )
      MODIFIED(
      SHORT(
      const BigReal* const slow_i = slow_table + 4*table_i;
      slow_a +=    modf_mod * slow_i[0];
      slow_b += 2.*modf_mod * slow_i[1];
      slow_c += 4.*modf_mod * slow_i[2];
      slow_d += 6.*modf_mod * slow_i[3];
      )
      NOSHORT(
      slow_d -= modf_mod * table_four_i[12];
      slow_c -= modf_mod * table_four_i[13];
      slow_b -= modf_mod * table_four_i[14];
      slow_a -= modf_mod * table_four_i[15];
      )
      )
      slow_d *= kqq;
      slow_c *= kqq;
      slow_b *= kqq;
      slow_a *= kqq;

      ENERGY(
      register BigReal slow_val =
	( ( diffa * slow_d *(1/6.)+ slow_c * (1/4.)) * diffa + slow_b *(1/2.)) * diffa + slow_a;
      fullElectEnergy -= LAM(lambda_pair *) slow_val;
      FEP( fullElectEnergy_s -= d_lambda_pair * slow_val; )
      )

      INT( {
      register BigReal slow_dir =
	( diffa * slow_d + slow_c ) * diffa + slow_b;
      reduction[pairElectForceIndex_X] += force_sign * slow_dir * p_ij_x;
      reduction[pairElectForceIndex_Y] += force_sign * slow_dir * p_ij_y;
      reduction[pairElectForceIndex_Z] += force_sign * slow_dir * p_ij_z;
      } )

      FAST(
      NOSHORT(
      slow_d += vdw_d;
      slow_c += vdw_c;
      slow_b += vdw_b;
      slow_a += vdw_a;
      )
      )

      register BigReal slow_dir =
	( diffa * slow_d + slow_c ) * diffa + slow_b;
      BigReal fullforce_r = slow_dir LAM(* lambda_pair);

      {
      register BigReal tmp_x = fullforce_r * p_ij_x;
      PAIR( fullElectVirial_xx += tmp_x * p_ij_x; )
      PAIR( fullElectVirial_xy += tmp_x * p_ij_y; )
      PAIR( fullElectVirial_xz += tmp_x * p_ij_z; )
      fullf_i_x += tmp_x;
      fullf_1[j].x -= tmp_x;
      register BigReal tmp_y = fullforce_r * p_ij_y;
      PAIR( fullElectVirial_yy += tmp_y * p_ij_y; )
      PAIR( fullElectVirial_yz += tmp_y * p_ij_z; )
      fullf_i_y += tmp_y;
      fullf_1[j].y -= tmp_y;
      register BigReal tmp_z = fullforce_r * p_ij_z;
      PAIR( fullElectVirial_zz += tmp_z * p_ij_z; )
      fullf_i_z += tmp_z;
      fullf_1[j].z -= tmp_z;

      PPROF(
        const BigReal p_j_z = p_j->position.z;
        int n2 = (int)floor((p_j_z-pressureProfileMin)*invThickness);
        pp_clamp(n2, pressureProfileSlabs);
        int p_j_partition = p_j->partition;

        pp_reduction(pressureProfileSlabs, n1, n2, 
                     p_i_partition, p_j_partition, pressureProfileAtomTypes,
                     tmp_x*p_ij_x, tmp_y * p_ij_y, tmp_z*p_ij_z,
                     pressureProfileReduction);

      )

      }
      )

   } // for pairlist

