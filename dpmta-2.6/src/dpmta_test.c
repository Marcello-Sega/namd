
/* 
*  dpmta_test.c - master PVM process for DPMTA implentation
*
*  w. t. rankin
*
*  Copyright (c) 1994 Duke University
*  All rights reserved
*
*  this program will instantiate and initialize the slave processes
*  and collect the resulting particle information for the parallel 
*  FMA algorithm implementation, utilizing the central master call
*  PMTAforce().
*
*  this program is an example of how to use the PMTA interface routines
*  instantiating and calling the PMTA slave processes from a central master
*  process.
*
*/

static char rcsid[]="$Id: dpmta_test.c,v 1.3 1997/09/29 23:58:44 jim Exp $";

/*
 * revision history:
 *
 * $Log: dpmta_test.c,v $
 * Revision 1.3  1997/09/29 23:58:44  jim
 * Incorporated changes from version 2.6.1 of DPMTA.
 *   - fixes for bad handling of empty/invalid multipoles when
 *     using large processor sets.
 *   - moved functions that provide data mapping to processors.  master
 *     and slave routines now call the same function in dpmta_distmisc.c
 * Also, switched pvmc.h back to pvm3.h.
 *
 * Revision 1.2  1997/09/12 22:56:35  jim
 * Modifications to work with converse pvm.
 *
 * Revision 1.1  1997/09/05 19:42:10  jim
 * Original distribution.
 *
 * Revision 2.18  1997/05/07  21:27:57  chumphre
 * cubectr changed to cellctr and other minor changes related to the interface
 *
 * Revision 2.17  1997/05/07  19:31:49  chumphre
 * implement a uniform interface for rectangular and ||-piped cells
 *
 * Revision 2.16  1997/05/07  18:59:46  chumphre
 * major changes to implement ||-piped cells
 *
 * Revision 2.15  1997/05/06  17:05:05  wrankin
 * - dpmta_test now generates a random boxsize and calls PMTAresize
 *   on every cycle.
 * - dpmta_test now has a -DSIZECHECK flag that prints ou the size of
 *   the serial DPMTA process every iteration (for debugging purposes
 *
 * Revision 2.14  1997/03/12  20:15:16  wrankin
 * updates to timer codes
 *
 * Revision 2.13  1997/02/26  16:54:34  wrankin
 * added support for serial non-pvm compilation of DPMTA
 *
 * Revision 2.12  1996/11/18  19:29:44  wrankin
 * update to virial code to use multipole virial calculation
 * old virial code is turned on by defining OLDVIRIAL
 * direct test program now computes absolute virial tensor.
 *
 * Revision 2.11  1996/10/29  19:34:16  wrankin
 * changes to makefile to support other platforms
 * fix for multi-master communication interface
 * new global dumpdata routine to support direct pbc/macro verification
 *
 * Revision 2.10  1996/10/28  23:02:01  wrankin
 * additions to test routines to provide processing paramenters as
 *   part of output file.
 * dpmta_direct will now perform macroscopic computations, reading
 *   in processing parameters from particle position file.
 *
 * Revision 2.9  1996/09/24  18:44:01  wrankin
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
 * Revision 2.8  1996/08/20  17:12:56  wrankin
 * Enhancement to allow any number of processes to be utilized.
 * Formally the number of precesses was limited to a power of 2.
 *
 * Revision 2.7  1996/08/14  16:07:01  wrankin
 * Fixed LJ computations
 *
 * Revision 2.6  1996/08/09  15:31:17  chumphre
 * Implemented noncubic (rectangular) boxes.
 *
 * Revision 2.5  1996/02/29  21:14:02  wrankin
 * New relaease: 2.4 (.1)
 *    - simplified calling structure for initialization
 *    - macroscopic periodic code
 *    - fixed PBC calculation - all particles are now stored in the
 *      cell as positions relative to the cell center.
 *    - virial preasure tensor computed
 *    - fix to allow particles on the outer cube boundary to be
 *      included. (UNC fix)
 *    - fix to order reception of particle data during the distributed
 *      calling sequence (UIUC fix)
 *    - removed M2L code that didn't use transfer matrices
 *    - early hooks in to perform interaction list sorting.
 *    - fixed LJ scaling factor for 1/r^12 potential.
 *    - cleaned up the LJ interface.
 *    - and of course, my continued efforts to ANSI-fy this beast.
 *
 * Revision 2.4  1996/02/12  15:21:15  wrankin
 * Added Macroscopic Assemblies code.
 *
 * Revision 2.3  1995/11/29  22:29:39  wrankin
 * addition of periodic boundary code
 * beginning work on virial code.
 *
 * Revision 2.2  1995/10/01  21:46:28  wrankin
 * implemented LJ multipole calculations - compile time option only
 *
 * Revision 2.1  1995/06/13  04:26:21  wrankin
 * No updates.
 * Sources placed under CVS.
 *
 * Revision 1.7  1995/02/24  21:23:32  wrankin
 * added interaction list parameter theta as a command line argument
 *
 * Revision 1.6  1995/01/13  17:05:58  wrankin
 * new particle distributions are generated on each iteration
 *
 * Revision 1.5  1994/12/06  19:15:49  wrankin
 * moved center of simulation cell to a non-zero value
 * changed size of simulation cell
 * cleaned up data dump output
 *
 * Revision 1.4  1994/11/30  18:07:35  wrankin
 * changed data structures to use "dpmta.h" types.
 *
 * Revision 1.3  1994/11/28  04:27:53  wrankin
 * added hooks to handle multiple slave calls to PMTAforce
 *
 * Revision 1.2  1994/11/17  00:27:52  wrankin
 * added code to handled variable cell center
 *
 * Revision 1.1  1994/10/26  02:51:55  wrankin
 * Initial revision
 *
 *
*/

/* include files */
#include <stdio.h>
#include <stdlib.h>
#ifndef SERIAL
#include "pvm3.h"
#endif
#include "dpmta.h"

#ifdef TIME
#include <sys/times.h>
#if defined HPPA
   #include<time.h>
#elif defined LINUX
   #include <sys/time.h>
#elif defined CRAY
   #include <sys/time.h>
#elif defined SGI64
   #include <sys/time.h>
#elif defined RS6K
   #include <sys/time.h>
   #define CLK_TCK 100
#endif
struct timeval trun;
double elapse_time;
#endif

#ifdef SERIAL
#ifdef SIZECHECK
#include <unistd.h>
#include <sys/time.h>
#include <sys/resource.h>
struct rusage sysdata;
int rusage_rtn;
#endif
#endif

main( int argc, char *argv[] ) {

   int i, j, k;                  /* loop counters */
   int mytid;                    /* master's task id */
   int iter;                     /* number of iterations */
   int tids[64];                 /* tids of pmta slaves */
   int num_parts;                /* number of parts to simulate */
   PmtaInitData initdata;        /* DPMTA initialization data */
   PmtaParticlePtr particlelist; /* array of particle information */
   PmtaPartInfoPtr forcelist;    /* array of resulting force information */
   PmtaPartInfoPtr flist_lj;     /* array of resulting LJ force info */
#ifdef PIPED
   double rn1, rn2, rn3;         /* Random numbers used for particle generation */
#endif
 
#if defined VIRIAL || defined OLDVIRIAL
   double vpot;
   PmtaVector vfor;
   double vpot_lj;
   PmtaVector vfor_lj;
#endif

   /* check and load command line arguments */
   if (argc != 9) {
      fprintf(stderr,
         "%s <#procs> <#lvls> <#parts> <fft> <mp> <theta> <iter> <pbc>\n",
         argv[0]);
      exit(-1);
      }

   initdata.nprocs = atoi(argv[1]);
   initdata.nlevels = atoi(argv[2]);
   num_parts = atoi(argv[3]);
   initdata.fft = atoi(argv[4]);
   initdata.mp = atoi(argv[5]);
   initdata.theta = atof(argv[6]);
   iter = atoi(argv[7]);

   initdata.pbc = atoi(argv[8]);
   initdata.kterm = 0;
   if (initdata.pbc > 0) {
      initdata.kterm = initdata.pbc - 1;
      initdata.pbc = 1;
   }


   /* other constants */
   initdata.mp_lj = 4;
   initdata.fftblock = 4;
#ifndef PIPED
   initdata.v1.x = 32.0;
   initdata.v2.y = 32.0;
   initdata.v3.z = 32.0;
#else
   initdata.v1.x = 64.0;
   initdata.v1.y = 0.0;
   initdata.v1.z = 0.0;
   initdata.v2.x = 32.0;
   initdata.v2.y = 32.0;
   initdata.v2.z = 0.0;
   initdata.v3.x = 0.0;
   initdata.v3.y = 0.0;
   initdata.v3.z = 32.0;
#endif
   initdata.cellctr.x = 16.0;
   initdata.cellctr.y = 16.0;
   initdata.cellctr.z = 16.0;

   initdata.calling_num = 0;
   initdata.calling_tids = (int *)NULL;

#ifndef SERIAL
   /* enroll master in pvm */
   mytid = pvm_mytid();
#endif

   /* initialize slave tasks */
   PMTAinit(&(initdata), tids);


   /****************************************************************
   *
   *  particle list processing -
   *
   *  create the particle lists and send each processor the list of 
   *  particles in each cell.
   *
   */

   /* allocate particle arrays */
   particlelist = (PmtaParticlePtr)malloc(num_parts * sizeof(PmtaParticle));
   forcelist = (PmtaPartInfoPtr)malloc(num_parts * sizeof(PmtaPartInfo));
   flist_lj = (PmtaPartInfoPtr)malloc(num_parts * sizeof(PmtaPartInfo));

   for (i=0; i<iter; i++) { 

#ifndef PIPED
      initdata.v1.x = 32.0 + (drand48() * 8.0);
      initdata.v2.y = 32.0;
      initdata.v3.z = 32.0;

      initdata.cellctr.x = 16.0;
      initdata.cellctr.y = 16.0;
      initdata.cellctr.z = 16.0;
#endif

      PMTAresize(&(initdata.v1), &(initdata.v2), &(initdata.v3), 
                 &(initdata.cellctr));
     
      /* create all particles */
      for (j=0; j<num_parts; j++) {
#ifndef PIPED
	 particlelist[j].p.x = (drand48() - 0.5) * initdata.v1.x +
	    initdata.cellctr.x;
	 particlelist[j].p.y = (drand48() - 0.5) * initdata.v2.y +
	    initdata.cellctr.y;
	 particlelist[j].p.z = (drand48() - 0.5) * initdata.v3.z +
	    initdata.cellctr.z;
#else
         rn1 = drand48() - 0.5;         
         rn2 = drand48() - 0.5;         
         rn3 = drand48() - 0.5;         
         particlelist[j].p.x = (rn1 * initdata.v1.x + rn2 * initdata.v2.x +
		rn3 * initdata.v3.x) + initdata.cellctr.x;
         particlelist[j].p.y = (rn1 * initdata.v1.y + rn2 * initdata.v2.y +
		rn3 * initdata.v3.y) + initdata.cellctr.y;
         particlelist[j].p.z = (rn1 * initdata.v1.z + rn2 * initdata.v2.z +
		rn3 * initdata.v3.z) + initdata.cellctr.z;
#endif

	 if ( (j%2) == 0 )
	    particlelist[j].q = 1.0;
	 else
	    particlelist[j].q = -1.0;


	 /* not absolutely needed */
	 particlelist[j].a = 1.0;
	 particlelist[j].b = 0.0;

      } /* for j */

#ifdef TIME
      gettimeofday(&trun,0);
      elapse_time = (double)trun.tv_sec +
	 ((double)trun.tv_usec/1000000);
#endif
     
      /*
       *  send particles out to slaves and wait for results
       */

      PMTAforce(num_parts, particlelist, forcelist, flist_lj);

#ifdef SERIAL
#ifdef SIZECHECK
      /*
       * report the memory usage
       */
      rusage_rtn = getrusage(RUSAGE_SELF,&sysdata);
      fprintf(stderr," iter=%d rtn=%d maxrss=%ld\n",
	      i, rusage_rtn, sysdata.ru_maxrss);
#endif
#endif

#if defined VIRIAL || defined OLDVIRIAL
      PMTAvirial( &vpot, &vfor, &vpot_lj, &vfor_lj );
      fprintf(stderr,"Virial = %lf (%lf,%lf,%lf)\n",
         vpot, vfor.x, vfor.y, vfor.z);
#ifdef COMP_LJ
      fprintf(stderr,"LJ Virial = %lf (%lf,%lf,%lf)\n",
         vpot_lj, vfor_lj.x, vfor_lj.y, vfor_lj.z);
#endif
#endif
   } /* for i */

#ifdef TIME
      gettimeofday(&trun,0);
      elapse_time = (double)trun.tv_sec +
	 ((double)trun.tv_usec/1000000) - elapse_time;
      fprintf(stderr,"Total Elapse Time = %lg sec.\n\n", elapse_time);
#endif


#ifdef DATADUMP
   Dump_Results(num_parts, &(initdata), particlelist, forcelist, flist_lj);
#endif


   /*
   *  clean up and kill off slaves
   */

   PMTAexit();

#ifndef SERIAL
   pvm_exit();
#endif
   exit(0);
}

