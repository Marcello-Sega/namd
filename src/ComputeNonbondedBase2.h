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

#ifndef ARCH_POWERPC
#pragma ivdep
#endif
  
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

#if ( SHORT( FAST( 1+ ) ) 0 ) 
      Force *f_j = f_1 + j;
#endif
	
#if ( FULL( 1+ ) 0 )
      Force *fullf_j = fullf_1 + j;
#endif

      //Power PC aliasing and alignment constraints
#ifdef ARCH_POWERPC
      __alignx(16, table_four_i);
      FAST (
      __alignx(16, lj_pars);
      )
      __alignx(16, p_j);
      
#if ( FULL( 1+ ) 0 )
#pragma disjoint (*table_four_i, *fullf_j)
#pragma disjoint (*p_j,          *fullf_j)
#if ( SHORT( FAST( 1+ ) ) 0 ) 
#pragma disjoint (*f_j    , *fullf_j)
#pragma disjoint (*fullf_j, *f_j)
#endif   //Short + fast
#endif   //Full

#if ( SHORT( FAST( 1+ ) ) 0 ) 
#pragma disjoint (*table_four_i, *f_j)
#pragma disjoint (*p_j,          *f_j)
#pragma disjoint (*lj_pars,      *f_j)
      __prefetch_by_load ((void *)&f_j->x);
#endif //Short + Fast

#endif   //ARCH_POWERPC

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
        const BigReal lambda_1 = lambda_table_i[2*p_j->partition];
        const BigReal lambda_2 = lambda_table_i[2*p_j->partition+1];
        const BigReal lambda_vdw_1 = lambda_vdw_table_i[2*p_j->partition];
        const BigReal lambda_vdw_2 = lambda_vdw_table_i[2*p_j->partition+1];
        const BigReal lambda_elec_1 = lambda_elec_table_i[2*p_j->partition];
        const BigReal lambda_elec_2 = lambda_elec_table_i[2*p_j->partition+1];
      )

      LES( BigReal lambda_pair = lambda_table_i[p_j->partition]; )

      register const BigReal p_ij_x = p_i_x - p_j->position.x;
      register const BigReal p_ij_y = p_i_y - p_j->position.y;
      register const BigReal p_ij_z = p_i_z - p_j->position.z;
      
#if ( FAST(1+) 0 )
      const BigReal A = scaling * lj_pars->A;
      const BigReal B = scaling * lj_pars->B;

      BigReal vdw_d = A * table_four_i[0] - B * table_four_i[2];
      BigReal vdw_c = A * table_four_i[1] - B * table_four_i[3];
      BigReal vdw_b = A * table_four_i[4] - B * table_four_i[6];
      BigReal vdw_a = A * table_four_i[5] - B * table_four_i[7];

      FEP (
        // Yes this could be made faster using lookup tables, but 
        // we're aiming for clarity here...
        // Writing the equations in terms of real physical quantities
        // makes it easier for others to later adapt the FEP functions
        const BigReal r2 = p_ij_x*p_ij_x + p_ij_y*p_ij_y + p_ij_z*p_ij_z;

        // The VdW parameters with the _1 and _2 suffix correspond to 
        // lambda1 and lambda2, respectively, and are used to construct
        // the scaled FEP VdW potential
        const BigReal r2_1 = r2 + lambda_shift_table_i[2*p_j->partition];
        const BigReal r6_1 = r2_1*r2_1*r2_1;

        const BigReal r2_2 = r2 + lambda_shift_table_i[2*p_j->partition+1];
        const BigReal r6_2 = r2_2*r2_2*r2_2;
      )   
          
      ENERGY(
      NOT_FEP (
      register BigReal vdw_val =
        ( ( diffa * vdw_d * (1/6.)+ vdw_c * (1/4.)) * diffa + vdw_b *(1/2.)) * diffa + vdw_a;
      vdwEnergy -= LAM(lambda_pair *) vdw_val;
      )
      FEP(
        // switching function (this is correct whether switching is active or not)
        const BigReal switchmul = r2 > switchdist2? \
                 switchfactor*(cutoff2 - r2)*(cutoff2 - r2)*(cutoff2 - 3.*switchdist2 + 2.*r2) \
                 : 1.;
        
        const BigReal fep_vdw_energy = A/(r6_1*r6_1) - B/r6_1; // needed later
           
        // modified FEP potential for vdW
        vdwEnergy   += lambda_vdw_1 * fep_vdw_energy * switchmul;
        vdwEnergy_s += lambda_vdw_2 * (A/(r6_2*r6_2) - B/r6_2) * switchmul;
      )
      ) // ENERGY
#endif // FAST

#if ( FAST(1+) 0 )
      INT( 
      register BigReal vdw_dir =
      ( diffa * vdw_d + vdw_c ) * diffa + vdw_b;
      //BigReal force_r =  LAM(lambda_pair *) vdw_dir;
      reduction[pairVDWForceIndex_X] += force_sign * vdw_dir * p_ij_x;
      reduction[pairVDWForceIndex_Y] += force_sign * vdw_dir * p_ij_y;
      reduction[pairVDWForceIndex_Z] += force_sign * vdw_dir * p_ij_z;
      )

#if ( SHORT(1+) 0 ) // Short-range electrostatics

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
      
      NOT_FEP (
        electEnergy -=  LAM(lambda_pair *) fast_val;
      )
      FEP(
        electEnergy   -= lambda_elec_1 * fast_val;
        electEnergy_s -= lambda_elec_2 * fast_val; 
      )
      ) //ENERGY

      INT(
      register BigReal fast_dir =
      ( diffa * fast_d + fast_c ) * diffa + fast_b;
      // force_r -= -1.0 * LAM(lambda_pair *) fast_dir;
      reduction[pairElectForceIndex_X] +=  force_sign * fast_dir * p_ij_x;
      reduction[pairElectForceIndex_Y] +=  force_sign * fast_dir * p_ij_y;
      reduction[pairElectForceIndex_Z] +=  force_sign * fast_dir * p_ij_z;
      )
      }

      // Combined short-range electrostatics and VdW force:
      NOT_FEP(
        fast_d += vdw_d;
        fast_c += vdw_c;
        fast_b += vdw_b;
        fast_a += vdw_a;  // not used!
        register BigReal fast_dir =
                    (diffa * fast_d + fast_c) * diffa + fast_b;
        BigReal force_r =  LAM(lambda_pair *) fast_dir;
      )
      FEP(
        // Note: switching has such a minor effect on the magnitude of the
        // force (~10ppm) that's it's not clear whether it's worth computing...
        // but we do it anyways!
        const BigReal switchmul2 = (r2 > switchdist2)? \
                 12.*switchfactor*(cutoff2 - r2)*(r2 - switchdist2) : 0.;
      
        // FEP force for Coulomb and vdW
        register BigReal fast_elect_dir = (diffa * fast_d + fast_c) * diffa + fast_b;
        const BigReal force_r = lambda_elec_1 * fast_elect_dir \
                          +  lambda_vdw_1 * (  \
                          + (12.*fep_vdw_energy + 6.*B/r6_1)/r2_1 * switchmul \
                          + fep_vdw_energy * switchmul2);
      )
          
      register BigReal tmp_x = force_r * p_ij_x;
      PAIR( virial_xx += tmp_x * p_ij_x; )
      PAIR( virial_xy += tmp_x * p_ij_y; )
      PAIR( virial_xz += tmp_x * p_ij_z; )

      f_i_x += tmp_x;
      f_j->x -= tmp_x;

      register BigReal tmp_y = force_r * p_ij_y;
      PAIR( virial_yy += tmp_y * p_ij_y; )
      PAIR( virial_yz += tmp_y * p_ij_z; )
      f_i_y += tmp_y;
      f_j->y -= tmp_y;
      
      register BigReal tmp_z = force_r * p_ij_z;
      PAIR( virial_zz += tmp_z * p_ij_z; )
      f_i_z += tmp_z;
      f_j->z -= tmp_z;

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

#if ( FULL (EXCLUDED( SHORT ( 1+ ) ) ) 0 ) 
      const BigReal* const slow_i = slow_table + 4*table_i;

#ifdef ARCH_POWERPC  //Alignment and aliasing constraints
      __alignx (16, slow_i);
#if ( SHORT( FAST( 1+ ) ) 0 ) 
#pragma disjoint (*slow_i, *f_j)
#endif
#pragma disjoint (*slow_i, *fullf_j)
#endif  //ARCH_POWERPC

#endif //FULL 


#if ( FULL (MODIFIED( SHORT ( 1+ ) ) ) 0 ) 
      const BigReal* const slow_i = slow_table + 4*table_i;

#ifdef ARCH_POWERPC //Alignment and aliasing constraints
      __alignx (16, slow_i);
#if ( SHORT( FAST( 1+ ) ) 0 ) 
#pragma disjoint (*slow_i, *f_j)
#endif
#pragma disjoint (*slow_i, *fullf_j)
#endif //ARCH_POWERPC

#endif //FULL
      
      FULL(
      BigReal slow_d = table_four_i[8 SHORT(+ 4)];
      BigReal slow_c = table_four_i[9 SHORT(+ 4)];
      BigReal slow_b = table_four_i[10 SHORT(+ 4)];
      BigReal slow_a = table_four_i[11 SHORT(+ 4)];
      EXCLUDED(
      SHORT(
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
      
      NOT_FEP (
        fullElectEnergy -= LAM(lambda_pair *) slow_val;
      )
          
      FEP(
        // PME is not scaled, so we use "lambda", not lambda_elec
        fullElectEnergy   -= lambda_1 * slow_val;
        fullElectEnergy_s -= lambda_2 * slow_val; 
      )
      ) // ENERGY

      INT( {
      register BigReal slow_dir =
	( diffa * slow_d + slow_c ) * diffa + slow_b;
      reduction[pairElectForceIndex_X] += force_sign * slow_dir * p_ij_x;
      reduction[pairElectForceIndex_Y] += force_sign * slow_dir * p_ij_y;
      reduction[pairElectForceIndex_Z] += force_sign * slow_dir * p_ij_z;
      } )


      NOT_FEP ( 
        FAST(
        NOSHORT(
        slow_d += vdw_d;
        slow_c += vdw_c;
        slow_b += vdw_b;
        slow_a += vdw_a; // unused!
        )
        )
        register BigReal slow_dir = (diffa * slow_d + slow_c) * diffa + slow_b;
        BigReal fullforce_r = slow_dir LAM(* lambda_pair);
      )
      FEP ( 
        // PME is not scaled, so we use "lambda", not lambda_elec
        register BigReal slow_dir = (diffa * slow_d + slow_c) * diffa + slow_b;
        BigReal fullforce_r = lambda_1 * slow_dir;
          
        FAST( NOSHORT(
          const BigReal switchmul2 = (r2 > switchdist2)? \
              12.*switchfactor*(cutoff2 - r2)*(r2 - switchdist2) : 0.;

          fullforce_r += lambda_vdw_1 * ( \
              (12.*fep_vdw_energy + 6.*B/r6_1)/r2_1 * switchmul \
               + fep_vdw_energy * switchmul2);
        ))
      )
          
      {
      register BigReal tmp_x = fullforce_r * p_ij_x;
      PAIR( fullElectVirial_xx += tmp_x * p_ij_x; )
      PAIR( fullElectVirial_xy += tmp_x * p_ij_y; )
      PAIR( fullElectVirial_xz += tmp_x * p_ij_z; )
      fullf_i_x += tmp_x;
      fullf_j->x -= tmp_x;
      register BigReal tmp_y = fullforce_r * p_ij_y;
      PAIR( fullElectVirial_yy += tmp_y * p_ij_y; )
      PAIR( fullElectVirial_yz += tmp_y * p_ij_z; )
      fullf_i_y += tmp_y;
      fullf_j->y -= tmp_y;
      register BigReal tmp_z = fullforce_r * p_ij_z;
      PAIR( fullElectVirial_zz += tmp_z * p_ij_z; )
      fullf_i_z += tmp_z;
      fullf_j->z -= tmp_z;

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

