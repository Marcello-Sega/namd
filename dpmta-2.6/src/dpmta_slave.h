/****************************************************************
*
*  dpmta_slave.h - global variables for dpmta pvm implementation
*
*  w. t. rankin
* 
*  these are the external references to the  global variables
*  used by the slave process during computation
*
*  Copyright (c) 1994 Duke University
*  All rights reserved
*
*  RCS info: $Id: dpmta_slave.h,v 1.1 1997/09/05 19:41:57 jim Exp $
*/

/*
 *  revision history:
 *  $Log: dpmta_slave.h,v $
 *  Revision 1.1  1997/09/05 19:41:57  jim
 *  Original distribution.
 *
 * Revision 2.21  1997/05/07  18:59:30  chumphre
 * major changes to implement ||-piped cells
 *
 * Revision 2.20  1996/11/18  19:29:30  wrankin
 * update to virial code to use multipole virial calculation
 * old virial code is turned on by defining OLDVIRIAL
 * direct test program now computes absolute virial tensor.
 *
 * Revision 2.19  1996/09/24  18:42:09  wrankin
 * many changes for support of version 2.5
 *  - non cubic cells
 *  - non power-of-2 processors
 *  - fixed the resize code
 *  - new virial interface
 *  - changes to macroscopic code (still not working)
 *  - code cleanup
 *  - new test program that reads PDB file
 *  - update code for T3D support
 *
 * Revision 2.18  1996/08/20  17:12:38  wrankin
 * Enhancement to allow any number of processes to be utilized.
 * Formally the number of precesses was limited to a power of 2.
 *
 * Revision 2.17  1996/08/09  15:30:49  chumphre
 * Implemented noncubic (rectangular) boxes.
 *
 * Revision 2.16  1996/02/12  15:20:59  wrankin
 * Added Macroscopic Assemblies code.
 *
 * Revision 2.15  1995/11/29  22:29:17  wrankin
 * addition of periodic boundary code
 * beginning work on virial code.
 *
 * Revision 2.14  1995/10/01  21:46:01  wrankin
 * implemented LJ multipole calculations - compile time option only
 *
 * Revision 2.13  1995/09/18  19:02:38  wrankin
 * updates to move the scaling of the simulation box to the slave
 * process.
 *
 * Revision 2.12  1995/07/01  03:26:57  wrankin
 * initial cut at precomputation of mpe transfer functions.
 * this works for vanilla M2L.  it does not work for FFT enhancements.
 *
 * Revision 2.11  1995/06/27  14:20:20  wrankin
 * Implementation of relative interaction list algorithm.
 * Interaction list generation code removed from master library.
 * Cell Table restructured as a 2-D table indexed by level/cell
 * All direct computations are now performed in a single pass.
 * General cleanup and removal of old files.
 *
 * Revision 2.10  1995/05/16  04:02:16  wrankin
 * added code to allow resizing the sim cube on the fly
 * the calling process now scales the input to a unit cube
 * and scales the results back to the proper value
 * the slave assumes all data sent to it is in a unit cube
 * with the bottom corner at <0,0,0>
 *
 * Revision 2.9  1995/04/09  22:21:06  wrankin
 * renaming of all global variables to prevent conflicts for
 * port of dpmta library to the T3D.
 *
 * Revision 2.8  1995/03/06  07:45:08  wrankin
 * master now sends Theta parameter to all slaves
 *
 * Revision 2.7  1995/01/13  17:07:01  wrankin
 * removed global local expansion buffer
 * changed type of global mpe buffer
 *
 * Revision 2.6  1994/11/30  18:09:30  wrankin
 * added globals for distributed calling sequence
 *
 * Revision 2.5  1994/11/28  04:27:53  wrankin
 * added hooks to handle multiple slave calls to PMTAforce
 *
 * Revision 2.4  1994/11/17  00:27:52  wrankin
 * added code to handled variable cell center
 *
 * Revision 2.3  1994/10/14  04:48:13  wrankin
 * added fft blocking global variable
 * added duke copyright notice
 * removed debugging ifdefs
 *
 * Revision 2.2  1994/07/25  23:08:40  wrankin
 * fixed message passing in upward pass portion of code
 * fixed mpe messages that were only sending half the complex data
 *
 * Revision 2.1  1994/06/02  12:52:27  wrankin
 * No change.  Update version number for Release 2
 *
 * Revision 1.5  1994/03/31  13:32:58  wrankin
 * added new global CubeLength, passed in from master program
 *
 * Revision 1.4  1994/03/26  16:29:23  wrankin
 * added FFT and MP globals
 * added LclSize global
 *
 * Revision 1.3  1994/03/18  17:50:36  wrankin
 * added global debugging outptu file pointer
 *
 * Revision 1.2  1994/03/14  15:32:14  wrankin
 * change global cell table name to CellTbl
 * misc cleanup
 *
 * Revision 1.1  1994/03/11  19:47:57  wrankin
 * Initial revision
 *
 *
*/


extern CellPtrPtrPtr Dpmta_CellTbl;    /* pointer to cell table list */
extern IlistPtr      Dpmta_Intlist;    /* interaction list table */
extern HlistPtr      Dpmta_Hlist;      /* mpe transfer matrix table */
extern IIdata        Dpmta_IIlist;     /* inverse interaction data table */

/* PVM processing parameters */

extern int Dpmta_Pid;                  /* process id of slave */
extern int Dpmta_Nproc;                /* total number of processors */
extern int Dpmta_MyTid;                /* my task id */
extern int Dpmta_Tids[];               /* task ids of other slave processes */
extern int Dpmta_MasterTid;            /* task id of master process */

extern int Dpmta_CallingNum;           /* interface for slave particles */
extern int Dpmta_CallingTids[];        /* task ids of slave procs */

/* MPE processing parameters */

extern int Dpmta_NumLevels;            /* number of levels in spatial decomp */
extern int Dpmta_DownPassStart;        /* level where downward pass starts */
extern int Dpmta_FFT;                  /* FFT processing flag */
extern int Dpmta_PBC;                  /* periodic boundary cond. flag */
extern int Dpmta_Mp;                   /* # terms in the multipole exp (p) */
extern int Dpmta_FftBlock;             /* FFT Blocking size (b) */
extern int Dpmta_MpeSize;              /* # Complex in the multipole exp */
extern int Dpmta_LclSize;              /* # Complex in the local exp */

extern double Dpmta_Theta;             /* multipole acceptance parameter */
#ifndef PIPED
extern Vector Dpmta_CellLength;        /* length of simulation cube */
#else
extern Vector Dpmta_PVector1;          /* ||-piped vectors and magnitudes */
extern Vector Dpmta_PVector2;
extern Vector Dpmta_PVector3;
extern double Dpmta_PV1Mag;
extern double Dpmta_PV2Mag;
extern double Dpmta_PV3Mag;
#endif
extern double Dpmta_MaxCellLen;	       /* length of longest side */
extern Vector Dpmta_CellCenter;        /* position of simulation cube center */
extern int Dpmta_Resize;               /* flag to initiate cell resizing */

#ifdef COMP_LJ
extern int Dpmta_Mp_LJ;                /* # terms in mpe for LJ potential */
extern int Dpmta_MpeSize_LJ;           /* # Complex in the LJ multipole exp */
extern int Dpmta_LclSize_LJ;           /* # Complex in the LJ local exp */

extern MtypeLJ Dpmta_Temp_Mpe_LJ;      /* temporary LJ multipole exp buffer */
#endif

/* misc processing variables */

extern int Dpmta_Scell[];              /* index to first owned cell */
extern int Dpmta_Ecell[];              /* index to last owned cells */
extern int Dpmta_RMcell[];             /* index to number of mpe's xfer'd */
extern int Dpmta_RLcell[];             /* index to number of loc's xfer'd */

extern Mtype Dpmta_Temp_Mpe;           /* temporary multipole exp buffer */

#if defined VIRIAL || defined OLDVIRIAL
extern Real Dpmta_Vpot;                /* virital potential */
extern Vector Dpmta_Vf;                /* virital force summation */
#ifdef COMP_LJ
extern Real Dpmta_Vpot_LJ;             /* virital potential */
extern Vector Dpmta_Vf_LJ;             /* virital force summation */
#endif
#endif

extern int Dpmta_Power8[];             /* arrays to aid cell table indexing */
extern int Dpmta_LevelLocate[];        /* arrays to aid cell table indexing */

#ifdef MACROSCOPIC
extern int Dpmta_K;                    /* # of levels in macro expansion */
#endif
