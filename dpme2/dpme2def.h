/* 
 *  Abdulnour Toukmaji 
 *  Copyright (c) 1996,1997 Duke University
 *  All rights reserved
 */

#ifndef DPME2DEF_H
#define DPME2DEF_H

/* global defines */

#define PI  3.14159265358979323846
#define  SQRT_PI  1.7724538509055160273 /* mathematica 15 digits*/
#define TwoBySqrtPi  1.12837916709551 /* mathematica 15 digits*/

/* this defines the frequency of updating the Link_cell data */
/* application must experiment with what is an acceptable value */
#define UPDATE_TIME 2
/********************************************************/
/*  array max size will be replaced with dynamic alloc in future version */

/* max number of atoms a node owns */
#define npmax 7000
/* max number of neighbour atoms */
#define nomax 8000
/* max number of owned and nearby atoms */
#define namax  npmax+nomax
/* max number of neighbour atoms per owned atom */
#define nnmax  200
/* max number of owned bins */
#define nbmax 1000
/* max number of atoms to send to all neighbors */
#define nemax 9000
/* snd/rcv buf size , multiple of 5 */
#define nfmax 10000
/* max number of swaps per timestep */
#define nsmax 32

/* macro definitions */
/********************************************************/
/* are now folded in dpme2.h */

/********************************************************/
/*         NEW STRUCTURES                               */
typedef struct processor {
   int node;
   int nprocs; /* total num of processors */
   int ncube;
   }  MYPROC;
typedef struct newpmeparticle {
  double x;
  double y;
  double z;
  double cg;
  double  id; /* global number of atom, should be int but this mk easy */
}  Pme2Particle;
typedef struct nfftgrid { /* nfft grid pts in xyz */
  int x;
  int y;
  int z;
} Grid;
typedef struct boxvector { /* only to specify box */
  double x;
  double y;
  double z;
  double origin;
}  PmeBVector;
/********************************************************/
/*            OLD_Data_Structures (DPME)                */
/*******************************************************/
typedef struct pmevector {
   double x;
   double y;
   double z;
   }  PmeVector;
typedef struct { double r, i; } doublecomplex; /* used for fft routines */

#endif
