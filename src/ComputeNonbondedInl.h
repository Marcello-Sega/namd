/**
***  Copyright (c) 1995, 1996, 1997, 1998, 1999, 2000 by
***  The Board of Trustees of the University of Illinois.
***  All rights reserved.
**/

/*
   Common operations for ComputeNonbonded classes
*/

#ifndef COMPUTENONBONDEDINL_H
#define COMPUTENONBONDEDINL_H

#ifndef NAMD_RESTRICT
#define restrict
#endif

#include "ComputeNonbondedUtil.h"
#include "SimParameters.h"
#include "Node.h"
#include "Molecule.h"
#include "LJTable.h"
#include "ReductionMgr.h"
#include "ReserveArray.h"

#include "PressureProfile.h"

inline int pairlist_from_pairlist(BigReal cutoff2,
				  BigReal p_i_x, BigReal p_i_y, BigReal p_i_z,
				  const CompAtom *p_j,
				  const plint *list, int list_size, plint *newlist,
				  BigReal r2_delta, BigReal *r2list) {
  
  BigReal cutoff2_delta = cutoff2 + r2_delta;
  plint *nli = newlist;
  BigReal *r2i = r2list;

  if ( list_size <= 0) return 0;

  int g = 0;

#ifndef SIMPLE_PAIRLIST
  //***************************************************************
  //* 4-way unrolled and software-pipelined 
  //***************************************************************

  int jout = 0;
  if ( list_size > 16) {
    // prefetch
    int jcur0 = list[g];
    int jcur1 = list[g + 1];
    int jcur2 = list[g + 2];
    int jcur3 = list[g + 3];
    
    int j0, j1, j2, j3;
    
    register  BigReal pj_x_0, pj_x_1, pj_x_2, pj_x_3; 
    register  BigReal pj_y_0, pj_y_1, pj_y_2, pj_y_3; 
    register  BigReal pj_z_0, pj_z_1, pj_z_2, pj_z_3; 
    
    register BigReal t_0, t_1, t_2, t_3, r2_0, r2_1, r2_2, r2_3;
    
    pj_x_0 = p_j[jcur0].position.x;
    pj_x_1 = p_j[jcur1].position.x;  
    pj_x_2 = p_j[jcur2].position.x;  
    pj_x_3 = p_j[jcur3].position.x;  
    pj_y_0 = p_j[jcur0].position.y; 
    pj_y_1 = p_j[jcur1].position.y;  
    pj_y_2 = p_j[jcur2].position.y;  
    pj_y_3 = p_j[jcur3].position.y;  
    pj_z_0 = p_j[jcur0].position.z; 
    pj_z_1 = p_j[jcur1].position.z;
    pj_z_2 = p_j[jcur2].position.z; 
    pj_z_3 = p_j[jcur3].position.z;
    
    for ( g = 4 ; g < list_size - 4; g += 4 ) {
      // compute 1d distance, 4-way parallel
      
      //Save the previous iterations values, gives more flexibility 
      //to the compiler to schedule the loads and the computation
      j0   =   jcur0;           j1   =   jcur1;
      j2   =   jcur2;           j3   =   jcur3;
      
      jcur0  =  list[g    ];    jcur1  =  list[g + 1];
      jcur2  =  list[g + 2];    jcur3  =  list[g + 3];

#if defined(ARCH_POWERPC) & !defined(MEM_OPT_VERSION)
      __dcbt ((void *) &p_j[jcur0]);
#endif      

      //Compute X distance
      t_0   =  p_i_x - pj_x_0;   t_1   =  p_i_x - pj_x_1;
      t_2   =  p_i_x - pj_x_2;   t_3   =  p_i_x - pj_x_3;
      
      r2_0  =  t_0 * t_0 + r2_delta; 
      r2_1  =  t_1 * t_1 + r2_delta;
      r2_2  =  t_2 * t_2 + r2_delta;
      r2_3  =  t_3 * t_3 + r2_delta;
      
      //Compute y distance
      t_0    =  p_i_y - pj_y_0;   t_1    =  p_i_y - pj_y_1;
      t_2    =  p_i_y - pj_y_2;   t_3    =  p_i_y - pj_y_3;
      r2_0  +=  t_0 * t_0;        r2_1  +=  t_1 * t_1;
      r2_2  +=  t_2 * t_2;        r2_3  +=  t_3 * t_3;
      
      //compute z distance
      t_0    =  p_i_z - pj_z_0;   t_1    =  p_i_z - pj_z_1;
      t_2    =  p_i_z - pj_z_2;   t_3    =  p_i_z - pj_z_3;
      r2_0  +=  t_0 * t_0;        r2_1  +=  t_1 * t_1;
      r2_2  +=  t_2 * t_2;        r2_3  +=  t_3 * t_3;
      
      pj_x_0 = p_j[jcur0].position.x;
      pj_x_1 = p_j[jcur1].position.x;  
      pj_x_2 = p_j[jcur2].position.x;  
      pj_x_3 = p_j[jcur3].position.x;  
      pj_y_0 = p_j[jcur0].position.y; 
      pj_y_1 = p_j[jcur1].position.y;  
      pj_y_2 = p_j[jcur2].position.y;  
      pj_y_3 = p_j[jcur3].position.y;  
      pj_z_0 = p_j[jcur0].position.z; 
      pj_z_1 = p_j[jcur1].position.z;
      pj_z_2 = p_j[jcur2].position.z; 
      pj_z_3 = p_j[jcur3].position.z;
      
      bool test0, test1, test2, test3;
      
      test0 = ( r2_0   <   cutoff2_delta );
      test1 = ( r2_1   <   cutoff2_delta );
      test2 = ( r2_2   <   cutoff2_delta );
      test3 = ( r2_3   <   cutoff2_delta );
      
      int jout0, jout1, jout2, jout3;

      jout0 = jout;
      nli[ jout0 ]  = j0;         r2i[ jout0 ] = r2_0;
      jout += test0;              jout1 = jout;
      nli[ jout1 ]  = j1;         r2i[ jout1 ] = r2_1;
      jout += test1;              jout2 = jout;
      nli[ jout2 ]  = j2;         r2i[ jout2 ] = r2_2;
      jout += test2;              jout3 = jout;
      nli[ jout3 ]  = j3;         r2i[ jout3 ] = r2_3;

      jout += test3;
    }
    g -= 4;
  }

  nli += jout;
  r2i += jout;  
#endif

  int j2 = list[g];
  BigReal p_j_x = p_j[j2].position.x;
  BigReal p_j_y = p_j[j2].position.y;
  BigReal p_j_z = p_j[j2].position.z;
  while ( g < list_size ) {
    int j = j2;
    j2 = list[++g];
    BigReal t2 = p_i_x - p_j_x;
    BigReal r2 = t2 * t2 + r2_delta;
    p_j_x = p_j[j2].position.x;
    t2 = p_i_y - p_j_y;
    r2 += t2 * t2;
    p_j_y = p_j[j2].position.y;
    t2 = p_i_z - p_j_z;
    r2 += t2 * t2;
    p_j_z = p_j[j2].position.z;
    if ( r2 <= cutoff2_delta ) {
      *nli= j; ++nli;
      *r2i = r2; ++r2i;
    }
  }

  return nli - newlist;
}

// clear all
// define interaction type (pair or self)
#define NBPAIR	1
#define NBSELF	2

#endif // COMPUTENONBONDEDINL_H

