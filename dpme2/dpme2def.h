/* 
 *  Abdulnour Toukmaji 
 *  Copyright (c) 1996,1997 Duke University
 *  All rights reserved
 */


/* global defines */

#ifndef PI
#undef PI
#define PI  3.14159265358979323846
#endif
#define  SQRT_PI  1.7724538509055160273 /* mathematica 15 digits*/
#define TwoBySqrtPi  1.12837916709551 /* mathematica 15 digits*/

/* this defines the frequency of updating the Link_cell data */
/* application must experiment with what is an acceptable value */
#define UPDATE_TIME 12

/* this defines if how we alloc recipF, either w/numatoms or as nlocal */
#define SAVE_RCPF 1 /* 1 =nlocal+MEM_INC, 0= npmax */
#define SAVE_DIRF 1 /* 1=nlocal +MEM_INC, 0= npmax  */
#define SAVE_ADJF 1 /* 1=nlocal +MEM_INC, 0= npmax */
#define NO_NINT 1 /* if =1 use new way to get around doing NINT min_img (org=0) */
/********************************************************/
#define MEM_INC 100 /* how much extra to realloc to  directF,recipF,adjustF */
#define MAXPROC         16   /* the max number of slave processors (exclude master) */
/********************************************************/
/*  array max size will be replaced with dynamic alloc in future version */

/* max number of atoms a node owns */
#define npmax 5000

/* max number of neighbour atoms */
#define nomax 17000

/* max number of owned and nearby atoms */
#define namax  npmax+nomax

/* max number of neighbour atoms per owned atom */
#define nnmax  160

/* max number of owned bins */
#define nbmax 1000

/* max number of atoms to send to all neighbors */
#define nemax 15500

/* snd/rcv buf size , multiple of 5 */
#define nfmax 16500

/* max number of swaps per timestep */
#define nsmax 32

/* total nbrs for all atoms (nour added 4/97) */
#define npnmax npmax*nnmax

#ifdef __cplusplus
extern "C" {
#endif

/********************************************************/
/* macro definitions */
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
/********************************************************/
#define DPME_DEFINES_PMEVECTOR
typedef struct pmevector {
   double x;
   double y;
   double z;
   }  PmeVector;
typedef struct { double r, i; } doublecomplex; /* used for fft routines */

/********************************************************/
/*          Added to support the Master/Worker DPME2    */
/********************************************************/
typedef struct PeAtomPnt { /* processor point to atom list */
  int pnt;  /* pionter 1st atom this PE has */
  int tot; /* total atoms this PE has */
} PePnt;
/*--------------------------------------------------------------------------*/
typedef struct AtomInfoS {
  int numatoms;
  int atompnt, freepnt; /* index to 1st freespace , atom in list */
  int nlocal;/* num of atoms in a processor */
  int nother;/* number of nearby atoms that My node stores */
} AtomInfo;
/*--------------------------------------------------------------------------*/
typedef struct BoxInfoS {
  PmeBVector box; /* unit cell dimensions in x y z, origin */
  PmeVector prd;
  PmeVector prd2; /* contains xprd,yprd,zprd = box dimensions(half) */
  double cutoff; /* should be <= rs < Box/2 */
  double rs; /* this is the shell cutoff in multi-time step */
  PmeVector mc2;  /* box min. inner cutoff */
  double ewaldcof; /* ewald coefficient */
  double dtol; /* user defined tolerance parameter that affects ewaldcof value */
  double recip[9]; /* reciprocal box dimensions */
} BoxInfo;
/*--------------------------------------------------------------------------*/
typedef struct PeInfoS{
  MYPROC myproc; /* info about parallel config */
  int inst_node[MAXPROC]; /* tid # of PEs */
  int igrid; /* flag indicates if user will specify processor grid */
} PeInfo;
/*--------------------------------------------------------------------------*/
typedef struct BoundaryInfoS{
  int npdim[3]; /* # of processors in each dim */
  int need[3]; /*# boxes away in each dimension for swapping */
  double border[6]; /*boundaries of processors box in each Dim */
  int mpart[6]; /* node# of nbr proc in each dim */
}BndryInfo;
/*--------------------------------------------------------------------------*/
typedef struct BinInfoS{
  PmeVector nbin, binsize; /*  bin size in  each dimension */
  PmeVector mbin, mbinlo; /* # bins in my box and 1st bin addrs */
  double boundlo[nsmax]; /* atom-position boundaries */
  double boundhi [nsmax]; /* atom-position boundaries */
  int *bin, *binpnt; /* bins of atoms in Link_cell method */
} BinInfo;
/*--------------------------------------------------------------------------*/
typedef struct SwapInfoS{
   int slist[nemax]; /* the swap list of atoms to send out */
   int nslist[nsmax]; /* pointer to beginnig of swap list per swap */
   int *nlist; /* nbr list of my atoms=[npnmax] ,alloced dynamically */
   int nnlist[npmax]; /* pointers to start of nbr list of my atoms */
   int nswap; /* number of swaps in each direction */
   int spart[nsmax], rpart[nsmax]; /* node's to send/rcv data from */
   short int *imglist; /* min image vector of same size as nlist, alloc'd dynamic*/
}SwapInfo;
/*--------------------------------------------------------------------------*/
typedef struct GridInfoS{
  int order; /* order of the interpolation */
  Grid nfftgrd; /* the grid dimensions */
  int nfft; /* grid dim if cubic */
  double volume; /* box volume */
}GridInfo;  
/*************************************************************************/

#ifdef __cplusplus
}
#endif

