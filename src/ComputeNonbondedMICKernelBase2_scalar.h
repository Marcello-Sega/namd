#ifdef NAMD_MIC

  // For each entry in the pairlist...

  // DMK - DEBUG
  #define USE_POINTER_MATH_FOR_TABLES (0)

  // Auto-vectorize via pairlist padding
  #if __MIC_PAD_PLGEN_CTRL != 0

    #if MIC_HANDCODE_FORCE_SINGLE != 0
      const int _plI_fs_outer_step = 16;
    #else
      const int _plI_fs_outer_step = 8;
    #endif

    #pragma novector
    for (int _plI_fs_outer = 0; _plI_fs_outer < plSize; _plI_fs_outer += _plI_fs_outer_step) {

      // Preload i value here (use broadcast)...
      const int i = (plArray[_plI_fs_outer] >> 16) & 0xFFFF;

      // Preload x,y,z,q values here
      #if MIC_HANDCODE_FORCE_SOA_VS_AOS != 0
        const CALC_TYPE p_i_x = ((CALC_TYPE)p_0[i].x) + ((CALC_TYPE)params.offset.x);
        const CALC_TYPE p_i_y = ((CALC_TYPE)p_0[i].y) + ((CALC_TYPE)params.offset.y);
        const CALC_TYPE p_i_z = ((CALC_TYPE)p_0[i].z) + ((CALC_TYPE)params.offset.z);
        const CALC_TYPE p_i_q = (CALC_TYPE)(p_0[i].charge);
        const int p_i_vdwType = pExt_0[i].vdw_type;
      #else
        const CALC_TYPE p_i_x = ((CALC_TYPE)p_0_x[i]) + ((CALC_TYPE)params.offset.x);
        const CALC_TYPE p_i_y = ((CALC_TYPE)p_0_y[i]) + ((CALC_TYPE)params.offset.y);
        const CALC_TYPE p_i_z = ((CALC_TYPE)p_0_z[i]) + ((CALC_TYPE)params.offset.z);
        const CALC_TYPE p_i_q = (CALC_TYPE)(p_0_q[i]);
        const int p_i_vdwType = pExt_0_vdwType[i];
      #endif

      #if 0
        FAST(SHORT( double tmp_x_i_sum = 0.0; ))
        FAST(SHORT( double tmp_y_i_sum = 0.0; ))
        FAST(SHORT( double tmp_z_i_sum = 0.0; ))
        #if MIC_EXCL_CHECKSUM != 0
          double tmp_w_i_sum = 0.0;
        #endif
        FULL( double fulltmp_x_i_sum = 0.0; )
        FULL( double fulltmp_y_i_sum = 0.0; )
        FULL( double fulltmp_z_i_sum = 0.0; )
        #if 1
          #pragma novector
        #else
          #pragma ivdep
        #endif
      #else
        double tmp_x_i_sum = 0.0;
        double tmp_y_i_sum = 0.0;
        double tmp_z_i_sum = 0.0;
        double tmp_w_i_sum = 0.0;
        double fulltmp_x_i_sum = 0.0;
        double fulltmp_y_i_sum = 0.0;
        double fulltmp_z_i_sum = 0.0;
        #pragma simd reduction(+:tmp_x_i_sum, tmp_y_i_sum, tmp_z_i_sum, tmp_w_i_sum, \
                                 fulltmp_x_i_sum, fulltmp_y_i_sum, fulltmp_z_i_sum)
      #endif
      for (int _plI_fs_inner = 0; _plI_fs_inner < _plI_fs_outer_step; _plI_fs_inner++) {
        const int plI = _plI_fs_outer + _plI_fs_inner;
        //if (/*plI < plSize &&*/ plArray[plI] >= 0) {
        if ((plArray[plI] & 0xFFFF) != 0xFFFF) {

          //// DMK - DEBUG
          //printf("[MIC] :: %d, %s-%d/%d :: i:%d, j:%d\n",
          //       params.ppI,
          //       NORMAL("N") MODIFIED("M") EXCLUDED("E"), plI, plSize,
          //       (plArray[plI] >> 16) & 0xFFFF, plArray[plI] & 0xFFFF
          //      );
          //fflush(NULL);

  // Scalar version of the code
  #else

    // DMK - NOTE : These loop_count values are loose, lower-bound guesses on my part (TODO : verify & refine)
    #if (0 PAIR(+1))
      #pragma loop_count (1000)
    #elif (0 SELF(+1))
      #pragma loop_count (10000)
    #endif
    for (int plI = 0; plI < plSize; plI++) {  

  #endif


    // Load the particle indicies
    const int ij = plArray[plI];
    //#if __MIC_PAD_PLGEN_CTRL != 0
    //  if (ij == -1) { continue; }
    //#endif
    #if __MIC_PAD_PLGEN_CTRL != 0
      // NOTE: moved before this loop, to the start of the _plI_fs_outer loop's body
    #else
      const int i = (ij >> 16) & 0xFFFF;   // If pairlist padding, i = "uniform"  (i, -1)
    #endif
    const int j = (ij      ) & 0xFFFF;   //                    , j = "unique" (

    // TODO | FIXME - Spread these out throughout the loop body (if possible) and
    //   change based on AoS versus SoA
    #if MIC_PREFETCH_DISTANCE > 0
      const int pfIJ = plArray[plI + MIC_PREFETCH_DISTANCE];
      const int pfI = (pfIJ >> 16) & 0xFFFF;
      const int pfJ = (pfIJ      ) & 0xFFFF;
      _mm_prefetch((char*)(p_0_x + pfI), MIC_PREFETCH_HINT);
      _mm_prefetch((char*)(p_0_y + pfI), MIC_PREFETCH_HINT);
      _mm_prefetch((char*)(p_0_z + pfI), MIC_PREFETCH_HINT);
      _mm_prefetch((char*)(p_0_q + pfI), MIC_PREFETCH_HINT);
      _mm_prefetch((char*)(f_0_x + pfI), MIC_PREFETCH_HINT);
      _mm_prefetch((char*)(f_0_y + pfI), MIC_PREFETCH_HINT);
      _mm_prefetch((char*)(f_0_z + pfI), MIC_PREFETCH_HINT);
      _mm_prefetch((char*)(p_1_x + pfJ), MIC_PREFETCH_HINT);
      _mm_prefetch((char*)(p_1_y + pfJ), MIC_PREFETCH_HINT);
      _mm_prefetch((char*)(p_1_z + pfJ), MIC_PREFETCH_HINT);
      _mm_prefetch((char*)(p_1_q + pfJ), MIC_PREFETCH_HINT);
      _mm_prefetch((char*)(f_1_x + pfJ), MIC_PREFETCH_HINT);
      _mm_prefetch((char*)(f_1_y + pfJ), MIC_PREFETCH_HINT);
      _mm_prefetch((char*)(f_1_z + pfJ), MIC_PREFETCH_HINT);
    #endif

    // Load atom information
    #if MIC_HANDCODE_FORCE_SOA_VS_AOS != 0
      #if __MIC_PAD_PLGEN_CTRL != 0
        // NOTE: moved before this loop, to the start of the _plI_fs_outer loop's body
      #else
        const CALC_TYPE p_i_x = ((CALC_TYPE)p_0[i].x) + ((CALC_TYPE)params.offset.x);
        const CALC_TYPE p_i_y = ((CALC_TYPE)p_0[i].y) + ((CALC_TYPE)params.offset.y);
        const CALC_TYPE p_i_z = ((CALC_TYPE)p_0[i].z) + ((CALC_TYPE)params.offset.z);
      #endif
      const CALC_TYPE p_j_x = (CALC_TYPE)(p_1[j].x);  // Neighboring gather to be optimized
      const CALC_TYPE p_j_y = (CALC_TYPE)(p_1[j].y);
      const CALC_TYPE p_j_z = (CALC_TYPE)(p_1[j].z);
    #else
      #if __MIC_PAD_PLGEN_CTRL != 0
        // NOTE: moved before this loop, to the start of the _plI_fs_outer loop's body
      #else
        const CALC_TYPE p_i_x = ((CALC_TYPE)p_0_x[i]) + ((CALC_TYPE)params.offset.x);
        const CALC_TYPE p_i_y = ((CALC_TYPE)p_0_y[i]) + ((CALC_TYPE)params.offset.y);
        const CALC_TYPE p_i_z = ((CALC_TYPE)p_0_z[i]) + ((CALC_TYPE)params.offset.z);
      #endif
      const CALC_TYPE p_j_x = (CALC_TYPE)(p_1_x[j]);
      const CALC_TYPE p_j_y = (CALC_TYPE)(p_1_y[j]);
      const CALC_TYPE p_j_z = (CALC_TYPE)(p_1_z[j]);
    #endif

    // Load position deltas and r2
    CALC_TYPE p_ij_x = p_i_x - p_j_x;
    CALC_TYPE p_ij_y = p_i_y - p_j_y;
    CALC_TYPE p_ij_z = p_i_z - p_j_z;

    #if REFINE_PAIRLISTS != 0
    CALC_TYPE r2 = (CALC_TYPE)(r2Array[plI]);
    #else
    CALC_TYPE r2 = (p_ij_x * p_ij_x) + (p_ij_y * p_ij_y) + (p_ij_z * p_ij_z) + r2_delta;
    if (r2 < cutoff2_delta) {
    #endif

      // Count this interaction as part of the exclChecksum
      #if MIC_EXCL_CHECKSUM != 0
        #if (0 MODIFIED(+1) EXCLUDED(+1))
          #if MIC_HANDCODE_FORCE_SOA_VS_AOS != 0
            #if __MIC_PAD_PLGEN_CTRL != 0
              tmp_w_i_sum += 1.0;
            #else
              f_0[i].w += 1.0;
            #endif
            f_1[j].w += 1.0;
          #else
            #if __MIC_PAD_PLGEN_CTRL != 0
              tmp_w_i_sum += 1.0;
            #else
              f_0_w[i] += 1.0;
            #endif
            f_1_w[j] += 1.0;
          #endif
        #endif
      #endif

      //#if 0
      //  const int table_i = (int)(((*((unsigned long long int *)(&r2))) >> 46)) + r2_delta_expc;
      //#elif 1
        #if MIC_HANDCODE_FORCE_SINGLE != 0
	  //const int table_i = ((int)((__intel_castf32_u32(r2)) >> 17)) + r2_delta_expc;
          const unsigned int table_i = ((int)((__intel_castf32_u32(r2)) >> 17)) + r2_delta_expc;
        #else
          //const int table_i = ((int)((__intel_castf64_u64(r2)) >> 46)) + r2_delta_expc;
          const unsigned int table_i = ((int)((__intel_castf64_u64(r2)) >> 46)) + r2_delta_expc;
        #endif
      //#elif 1
      //  union byte_order { double d; int i[2]; };
      //  byte_order r2Int;
      //  r2Int.d = r2;
      //  // DMK - TODO | FIXME : This is hard coded for little-endian at the moment !!
      //  const int table_i = (r2Int.i[1] >> 14) + r2_delta_expc;
      //#else
      //  const int table_i = ((reinterpret_cast<int>((float)r2)) >> 17) + r2_delta_expc;
      //#endif
      //__assume(table_i >= 0);

      //// DMK - DEBUG
      //printf("[MIC] ::   table_i:%d\n", table_i);
      //fflush(NULL);

      #if MIC_HANDCODE_FORCE_CALCR2TABLE != 0
        // From ComputeNonbondedUtil.C                    Simplified:
        //   r2_base = r2_delta * (1 << (i/64))             r2_base = r2_delta * (1 << (i/64))
        //   r2_del = r2_base / 64.0;                       r2_del = r2_base / 64.0;
        //   r2 = r2_base - r2_delta + r2_del * (i%64)      r2_table[i] = r2_base - r2_delta + r2_del * (i%64) + r2_delta;
        //   r2_table[i] = r2 + r2_delta;                               = r2_base + r2_del * (i%64)
        // NOTE: For i = 0, r2_table[0] = r2_delta + (r2_delta / 64) * 0 = r2_delta, so there no need
        //   to special case if table_i = 0 then r2_table[0] = r2_delta (see ComputeNonbondedUtil.C:606)
        CALC_TYPE r2_base = r2_delta * (1 << (table_i >> 6)); // avoid original divide (table_i / 64)
        CALC_TYPE r2_del = r2_base * ((CALC_TYPE)0.015625f);  // avoid original divide (r2_base / 64)
        //CALC_TYPE r2_table_i = r2_base + r2_del * (table_i % 64);  // NOTE: removing '+ r2_delta - r2_delta'
        CALC_TYPE r2_table_i = r2_base + r2_del * (table_i & 0x3F); //(table_i % 64);  // NOTE: removing '+ r2_delta - r2_delta'
      #else
        CALC_TYPE r2_table_i = r2_table[table_i];
      #endif
      CALC_TYPE diffa = r2 - r2_table_i;
      #if USE_POINTER_MATH_FOR_TABLES != 0
        const CALC_TYPE * const table_four_i = SHORT(table_short) NOSHORT(table_noshort) + (16 * table_i);
      #else
        const CALC_TYPE * const table_four_ptr = SHORT(table_short) NOSHORT(table_noshort);
        const int table_four_idx = 16 * table_i;
      #endif

      // NOTE : These charge values are already scaled by
      //   sqrt(COULOMB * scaling * dielectric_1).  See HomePatch.C.
      #if MIC_HANDCODE_FORCE_SOA_VS_AOS != 0
        #if __MIC_PAD_PLGEN_CTRL != 0
          // NOTE: moved before this loop, to the start of the _plI_fs_outer loop's body
        #else
          const CALC_TYPE p_i_q = (CALC_TYPE)(p_0[i].charge);
        #endif
        const CALC_TYPE p_j_q = (CALC_TYPE)(p_1[j].charge);
      #else
        #if __MIC_PAD_PLGEN_CTRL != 0
          // NOTE: moved before this loop, to the start of the _plI_fs_outer loop's body
        #else
          const CALC_TYPE p_i_q = (CALC_TYPE)(p_0_q[i]);
        #endif
        const CALC_TYPE p_j_q = (CALC_TYPE)(p_1_q[j]);
      #endif
      CALC_TYPE kqq = p_i_q * p_j_q;

      #if (0 FAST(+1))

        #if MIC_HANDCODE_FORCE_SOA_VS_AOS != 0
          #if __MIC_PAD_PLGEN_CTRL != 0
            // NOTE: moved before this loop, to the start of the _plI_fs_outer loop's body
          #else
            int p_i_vdwType = pExt_0[i].vdw_type;
          #endif
          int p_j_vdwType = pExt_1[j].vdw_type;
        #else
          #if __MIC_PAD_PLGEN_CTRL != 0
            // NOTE: moved before this loop, to the start of the _plI_fs_outer loop's body
          #else
            int p_i_vdwType = pExt_0_vdwType[i];
          #endif
          int p_j_vdwType = pExt_1_vdwType[j];
        #endif

        //// DMK - DEBUG
        //printf("[MIC] ::   p_i_vdwType:%d, p_j_vdwType:%d\n", p_i_vdwType, p_j_vdwType);
        //fflush(NULL);

        #if 0
          const CALC_TYPE * const lj_pars_base = lj_table_base_ptr
            + (4 * (p_i_vdwType * lj_table_dim + p_j_vdwType)) // 4 CALC_TYPEs per entry: 2 normal, 2 modified
            MODIFIED(+ 2);
          CALC_TYPE A = scaling * lj_pars_base[0];
          CALC_TYPE B = scaling * lj_pars_base[1];
        #else
          const int lj_pars_offset = (4 * (p_i_vdwType * lj_table_dim + p_j_vdwType)) MODIFIED(+ 2);
          CALC_TYPE A = scaling * lj_table_base_ptr[lj_pars_offset    ];
          CALC_TYPE B = scaling * lj_table_base_ptr[lj_pars_offset + 1];
        #endif

        // 16x16 AoS table lookup with transpose
        #if USE_POINTER_MATH_FOR_TABLES != 0
          CALC_TYPE vdw_d = A * table_four_i[0] - B * table_four_i[4];
          CALC_TYPE vdw_c = A * table_four_i[1] - B * table_four_i[5];
          CALC_TYPE vdw_b = A * table_four_i[2] - B * table_four_i[6];
          CALC_TYPE vdw_a = A * table_four_i[3] - B * table_four_i[7];
        #else
          CALC_TYPE vdw_d = A * table_four_ptr[table_four_idx + 0] - B * table_four_ptr[table_four_idx + 4];
          CALC_TYPE vdw_c = A * table_four_ptr[table_four_idx + 1] - B * table_four_ptr[table_four_idx + 5];
          CALC_TYPE vdw_b = A * table_four_ptr[table_four_idx + 2] - B * table_four_ptr[table_four_idx + 6];
          CALC_TYPE vdw_a = A * table_four_ptr[table_four_idx + 3] - B * table_four_ptr[table_four_idx + 7];
        #endif

        #if (0 ENERGY(+1))
          CALC_TYPE vdw_val = ((diffa * vdw_d * (1/6.0) + vdw_c * (1/4.0)) * diffa + vdw_b * (1/2.0)) * diffa + vdw_a;
          vdwEnergy -= vdw_val;
          // DMK - TODO | FIXME : Apply vdw_val to FEP(vdwEnergy_s)
        #endif

        #if (0 SHORT(+1))

          #if (0 NORMAL(+1))
            #if USE_POINTER_MATH_FOR_TABLES != 0
              CALC_TYPE fast_d = kqq * table_four_i[8];
              CALC_TYPE fast_c = kqq * table_four_i[9];
              CALC_TYPE fast_b = kqq * table_four_i[10];
              CALC_TYPE fast_a = kqq * table_four_i[11];
            #else
              CALC_TYPE fast_d = kqq * table_four_ptr[table_four_idx +  8];
              CALC_TYPE fast_c = kqq * table_four_ptr[table_four_idx +  9];
              CALC_TYPE fast_b = kqq * table_four_ptr[table_four_idx + 10];
              CALC_TYPE fast_a = kqq * table_four_ptr[table_four_idx + 11];
            #endif
          #endif
          #if (0 MODIFIED(+1))
            CALC_TYPE modfckqq = (1.0 - modf_mod) * kqq;
            #if USE_POINTER_MATH_FOR_TABLES != 0
              CALC_TYPE fast_d = modfckqq * table_four_i[8];
              CALC_TYPE fast_c = modfckqq * table_four_i[9];
              CALC_TYPE fast_b = modfckqq * table_four_i[10];
              CALC_TYPE fast_a = modfckqq * table_four_i[11];
            #else
              CALC_TYPE fast_d = modfckqq * table_four_ptr[table_four_idx +  8];
              CALC_TYPE fast_c = modfckqq * table_four_ptr[table_four_idx +  9];
              CALC_TYPE fast_b = modfckqq * table_four_ptr[table_four_idx + 10];
              CALC_TYPE fast_a = modfckqq * table_four_ptr[table_four_idx + 11];
            #endif
          #endif

          #if (0 ENERGY(+1))
            CALC_TYPE fast_val = ((diffa * fast_d * (1/6.0) + fast_c * (1/4.0)) * diffa + fast_b * (1/2.0)) * diffa + fast_a;
            #if (0 NOT_ALCHPAIR(+1))
              electEnergy -= fast_val;
              // DMK - TODO | FIXME : Apply fast_val to FEP(electEnergy_s)
            #endif
          #endif

          #if (0 NOT_ALCHPAIR(+1))
            fast_d += vdw_d;
            fast_c += vdw_c;
            fast_b += vdw_b;
            fast_a += vdw_a;
          #endif

          CALC_TYPE fast_dir = (fast_d * diffa + fast_c) * diffa + fast_b;
          CALC_TYPE force_r = fast_dir;

          CALC_TYPE tmp_x = force_r * p_ij_x;
          PAIR( virial_xx += tmp_x * p_ij_x; )
          PAIR( virial_xy += tmp_x * p_ij_y; )
          PAIR( virial_xz += tmp_x * p_ij_z; )
          #if MIC_HANDCODE_FORCE_SOA_VS_AOS != 0
            #if __MIC_PAD_PLGEN_CTRL != 0
              tmp_x_i_sum += tmp_x;
            #else
              f_0[i].x += tmp_x;
            #endif
            f_1[j].x -= tmp_x;
          #else
            #if __MIC_PAD_PLGEN_CTRL != 0
              tmp_x_i_sum += tmp_x;
            #else
              f_0_x[i] += tmp_x;
            #endif
            f_1_x[j] -= tmp_x;
          #endif

          CALC_TYPE tmp_y = force_r * p_ij_y;
          PAIR( virial_yy += tmp_y * p_ij_y; )
          PAIR( virial_yz += tmp_y * p_ij_z; )
          #if MIC_HANDCODE_FORCE_SOA_VS_AOS != 0
            #if __MIC_PAD_PLGEN_CTRL != 0
              tmp_y_i_sum += tmp_y;
            #else
              f_0[i].y += tmp_y;   /// Move out after inner loop
            #endif
            f_1[j].y -= tmp_y;
          #else
            #if __MIC_PAD_PLGEN_CTRL != 0
              tmp_y_i_sum += tmp_y;
            #else
              f_0_y[i] += tmp_y;
            #endif
            f_1_y[j] -= tmp_y;
          #endif

          CALC_TYPE tmp_z = force_r * p_ij_z;
          PAIR( virial_zz += tmp_z * p_ij_z; )
          #if MIC_HANDCODE_FORCE_SOA_VS_AOS != 0
            #if __MIC_PAD_PLGEN_CTRL != 0
              tmp_z_i_sum += tmp_z;
            #else
              f_0[i].z += tmp_z;
            #endif
            f_1[j].z -= tmp_z;
          #else
            #if __MIC_PAD_PLGEN_CTRL != 0
              tmp_z_i_sum += tmp_z;
            #else
              f_0_z[i] += tmp_z;
            #endif
            f_1_z[j] -= tmp_z;
          #endif

        #endif // SHORT
      #endif // FAST

      #if (0 FULL(+1))

        #if USE_POINTER_MATH_FOR_TABLES != 0
          CALC_TYPE slow_d = table_four_i[ 8 SHORT(+ 4)];
          CALC_TYPE slow_c = table_four_i[ 9 SHORT(+ 4)];
          CALC_TYPE slow_b = table_four_i[10 SHORT(+ 4)];
          CALC_TYPE slow_a = table_four_i[11 SHORT(+ 4)];
        #else
          CALC_TYPE slow_d = table_four_ptr[table_four_idx +  8 SHORT(+ 4)];
          CALC_TYPE slow_c = table_four_ptr[table_four_idx +  9 SHORT(+ 4)];
          CALC_TYPE slow_b = table_four_ptr[table_four_idx + 10 SHORT(+ 4)];
          CALC_TYPE slow_a = table_four_ptr[table_four_idx + 11 SHORT(+ 4)];
        #endif

        #if (0 SHORT( EXCLUDED(+1) MODIFIED(+1) ))
          #if USE_POINTER_MATH_FOR_TABLES != 0
            const CALC_TYPE * const slow_i = slow_table + 4 * table_i;
          #else
            const int slow_idx = 4 * table_i;
          #endif
        #endif
        #if (0 EXCLUDED(+1))
          #if (0 SHORT(+1))
            #if USE_POINTER_MATH_FOR_TABLES != 0
              slow_a += 1.0 * slow_i[3];  // AoS transpose (4 members)
              slow_b += 2.0 * slow_i[2];
              slow_c += 4.0 * slow_i[1];
              slow_d += 6.0 * slow_i[0];
            #else
              slow_a += 1.0 * slow_table[slow_idx + 3];  // AoS transpose (4 members)
              slow_b += 2.0 * slow_table[slow_idx + 2];
              slow_c += 4.0 * slow_table[slow_idx + 1];
              slow_d += 6.0 * slow_table[slow_idx + 0];
            #endif
          #endif
          #if (0 NOSHORT(+1))
            #if USE_POINTER_MATH_FOR_TABLES != 0
              slow_d -= table_four_i[12];
              slow_c -= table_four_i[13];
              slow_b -= table_four_i[14];
              slow_a -= table_four_i[15];
            #else
              slow_d -= table_four_ptr[table_four_idx + 12];
              slow_c -= table_four_ptr[table_four_idx + 13];
              slow_b -= table_four_ptr[table_four_idx + 14];
              slow_a -= table_four_ptr[table_four_idx + 15];
            #endif
          #endif
        #endif
        #if (0 MODIFIED(+1))
          #if (0 SHORT(+1))
            #if USE_POINTER_MATH_FOR_TABLES != 0
              slow_a += 1.0 * modf_mod * slow_i[3];
              slow_b += 2.0 * modf_mod * slow_i[2];
              slow_c += 4.0 * modf_mod * slow_i[1];
              slow_d += 6.0 * modf_mod * slow_i[0];
            #else
              slow_a += 1.0 * modf_mod * slow_table[slow_idx + 3];
              slow_b += 2.0 * modf_mod * slow_table[slow_idx + 2];
              slow_c += 4.0 * modf_mod * slow_table[slow_idx + 1];
              slow_d += 6.0 * modf_mod * slow_table[slow_idx + 0];
            #endif
          #endif
          #if (0 NOSHORT(+1))
            #if USE_POINTER_MATH_FOR_TABLES != 0
              slow_d -= modf_mod * table_four_i[12];
              slow_c -= modf_mod * table_four_i[13];
              slow_b -= modf_mod * table_four_i[14];
              slow_a -= modf_mod * table_four_i[15];
            #else
              slow_d -= modf_mod * table_four_ptr[table_four_idx + 12];
              slow_c -= modf_mod * table_four_ptr[table_four_idx + 13];
              slow_b -= modf_mod * table_four_ptr[table_four_idx + 14];
              slow_a -= modf_mod * table_four_ptr[table_four_idx + 15];
            #endif
          #endif
        #endif
        slow_d *= kqq;
        slow_c *= kqq;
        slow_b *= kqq;
        slow_a *= kqq;

        #if (0 ENERGY(+1))
          CALC_TYPE slow_val = ((diffa * slow_d * (1/6.0) + slow_c * (1/4.0)) * diffa + slow_b * (1/2.0)) * diffa + slow_a;
          #if (0 NOT_ALCHPAIR(+1))
            fullElectEnergy -= slow_val;
            // DMK - TODO | FIXME : Apply slow_val to FEP(fullElectEnergy_s)
          #endif
        #endif

        #if (0 NOT_ALCHPAIR(FAST(NOSHORT(+1))))
          slow_d += vdw_d;
          slow_c += vdw_c;
          slow_b += vdw_b;
          slow_a += vdw_a;
        #endif

        CALC_TYPE slow_dir = (diffa * slow_d + slow_c) * diffa + slow_b;
        CALC_TYPE fullforce_r = slow_dir;

        CALC_TYPE fulltmp_x = fullforce_r * p_ij_x;
        PAIR( fullElectVirial_xx += fulltmp_x * p_ij_x; )
        PAIR( fullElectVirial_xy += fulltmp_x * p_ij_y; )
        PAIR( fullElectVirial_xz += fulltmp_x * p_ij_z; )
        #if MIC_HANDCODE_FORCE_SOA_VS_AOS != 0
          #if __MIC_PAD_PLGEN_CTRL != 0
            fulltmp_x_i_sum += fulltmp_x;
          #else
            fullf_0[i].x += fulltmp_x;
          #endif
          fullf_1[j].x -= fulltmp_x;
        #else
          #if __MIC_PAD_PLGEN_CTRL != 0
            fulltmp_x_i_sum += fulltmp_x;
          #else
            fullf_0_x[i] += fulltmp_x;
          #endif
          fullf_1_x[j] -= fulltmp_x;
        #endif

        CALC_TYPE fulltmp_y = fullforce_r * p_ij_y;
        PAIR( fullElectVirial_yy += fulltmp_y * p_ij_y; )
        PAIR( fullElectVirial_yz += fulltmp_y * p_ij_z; )
        #if MIC_HANDCODE_FORCE_SOA_VS_AOS != 0
          #if __MIC_PAD_PLGEN_CTRL != 0
            fulltmp_y_i_sum += fulltmp_y;
          #else
            fullf_0[i].y += fulltmp_y;
          #endif
          fullf_1[j].y -= fulltmp_y;
        #else
          #if __MIC_PAD_PLGEN_CTRL != 0
            fulltmp_y_i_sum += fulltmp_y;
          #else
            fullf_0_y[i] += fulltmp_y;
          #endif
          fullf_1_y[j] -= fulltmp_y;
        #endif

        CALC_TYPE fulltmp_z = fullforce_r * p_ij_z;
        PAIR( fullElectVirial_zz += fulltmp_z * p_ij_z; )
        #if MIC_HANDCODE_FORCE_SOA_VS_AOS != 0
          #if __MIC_PAD_PLGEN_CTRL != 0
            fulltmp_z_i_sum += fulltmp_z;
          #else
            fullf_0[i].z += fulltmp_z;
          #endif
          fullf_1[j].z -= fulltmp_z;
        #else
          #if __MIC_PAD_PLGEN_CTRL != 0
            fulltmp_z_i_sum += fulltmp_z;
          #else
            fullf_0_z[i] += fulltmp_z;
          #endif
          fullf_1_z[j] -= fulltmp_z;
        #endif

      #endif // FULL

      // DMK - DEBUG
      #if (0 FAST(SHORT(+1))) && 0

      #if (0 FULL(+1))
	const char * const idStr = "10";
      #else
	const char * const idStr = "11";
      #endif

      if (params.p1 == 0 && params.p2 <= 1) {
	printf("[%s], %03d, %03d, %03d, %03d, 0, Interaction...\n",
               idStr, params.p1, params.p2, i, j
	      );
	printf("[%s], %03d, %03d, %03d, %03d, 1,   p_0:{%.8lf %.8lf %.8lf}, p_1:{%.8lf %.8lf %.8lf}\n",
               idStr, params.p1, params.p2, i, j,
               (p_0_x[i] + params.patch1_center_x), (p_0_y[i] + params.patch1_center_y), (p_0_z[i] + params.patch1_center_z),
               (p_1_x[j] + params.patch2_center_x), (p_1_y[j] + params.patch2_center_y), (p_1_z[j] + params.patch2_center_z)
              );
	printf("[%s], %03d, %03d, %03d, %03d, 1,   (raw) p_0:{%.8lf %.8lf %.8lf}, patch1_center:{%.8lf %.8lf %.8lf}\n",
               idStr, params.p1, params.p2, i, j,
               p_0_x[i], p_0_y[i], p_0_z[i],
               params.patch1_center_x, params.patch1_center_y, params.patch1_center_z
              );
	printf("[%s], %03d, %03d, %03d, %03d, 1,   (raw) p_1:{%.8lf %.8lf %.8lf}, patch2_center:{%.8lf %.8lf %.8lf}\n",
               idStr, params.p1, params.p2, i, j,
               p_1_x[j], p_1_y[j], p_1_z[j],
               params.patch2_center_x, params.patch2_center_y, params.patch2_center_z
              );
	printf("[%s], %03d, %03d, %03d, %03d, 2,   r2_delta:%.18le\n",
               idStr, params.p1, params.p2, i, j,
               r2_delta
              );
	printf("[%s], %03d, %03d, %03d, %03d, 3,   offset:{%.18lf %.18lf %.18lf}\n",
               idStr, params.p1, params.p2, i, j,
               params.offset.x, params.offset.y, params.offset.z
              );
	printf("[%s], %03d, %03d, %03d, %03d, 4,   p_ij:{%.18lf %.18lf %.18lf}\n",
               idStr, params.p1, params.p2, i, j,
               p_ij_x, p_ij_y, p_ij_z
              );
        printf("[%s], %03d, %03d, %03d, %03d, 5,   diffa:%.18le, fast_dcba:{ %.18le %.18le %.18le %.18le }\n",
               idStr, params.p1, params.p2, i, j,
               diffa, fast_d, fast_c, fast_b, fast_a
	      );
        printf("[%s], %03d, %03d, %03d, %03d, 5,   r2:%.18le, table_i:%d, force_r:%.18le\n",
               idStr, params.p1, params.p2, i, j,
               r2, table_i, force_r
              );
        printf("[%s], %03d, %03d, %03d, %03d, 5,   A:%.18le, B:%.18le\n",
               idStr, params.p1, params.p2, i, j,
               A, B
              );
	printf("[%s], %03d, %03d, %03d, %03d, 5,   tmp:{%.18lf %.18lf %.18lf}\n",
               idStr, params.p1, params.p2, i, j,
               tmp_x, tmp_y, tmp_z
              );
        printf("[%s], %03d, %03d, %03d, %03d, 6,   f_1:{%.18le %.18le %.18le}\n",
               idStr, params.p1, params.p2, i, j,
               f_1_x[j], f_1_y[j], f_1_z[j]
              );
      }
      #endif

    #if REFINE_PAIRLISTS == 0
    } // end if (r2 < cutoff2_delta)
    #endif

  // End of loops auto-vectorized via pairlist padding
  #if __MIC_PAD_PLGEN_CTRL != 0

      } // end if
    } // end for

    // Apply force reductions
    #if MIC_HANDCODE_FORCE_SOA_VS_AOS != 0
      #if MIC_EXCL_CHECKSUM != 0
        f_0[i].w += tmp_w_i_sum;
      #endif
    #else
      #if MIC_EXCL_CHECKSUM != 0
        f_0_w[i] += tmp_w_i_sum;
      #endif
    #endif

    #if (0 FAST(SHORT(+1)))
      #if MIC_HANDCODE_FORCE_SOA_VS_AOS != 0
        f_0[i].x += tmp_x_i_sum;
        f_0[i].y += tmp_y_i_sum;
        f_0[i].z += tmp_z_i_sum;
      #else
        f_0_x[i] += tmp_x_i_sum;
        f_0_y[i] += tmp_y_i_sum;
        f_0_z[i] += tmp_z_i_sum;
      #endif
    #endif

    #if (0 FULL(+1))
      #if MIC_HANDCODE_FORCE_SOA_VS_AOS != 0
        fullf_0[i].x += fulltmp_x_i_sum;
        fullf_0[i].y += fulltmp_y_i_sum;
        fullf_0[i].z += fulltmp_z_i_sum;
      #else
        fullf_0_x[i] += fulltmp_x_i_sum;
        fullf_0_y[i] += fulltmp_y_i_sum;
        fullf_0_z[i] += fulltmp_z_i_sum;
      #endif
    #endif

  } // end for


  // End of scalar loop
  #else

  } // end pairlist-loop

  #endif

  // DMK - DEBUG
  #undef USE_POINTER_MATH_FOR_TABLES

#endif  // NAMD_MIC
