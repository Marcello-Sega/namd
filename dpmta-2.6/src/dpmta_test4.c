
/* 
*  dpmta_test4.c - master PVM process for DPMTA implentation
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
*  the input for this program is a PDB file with a number of water molecules
*  in it.  the simulation region is assumed to be cubic.
*
*/

static char rcsid[]="$Id: dpmta_test4.c,v 1.1 1997/09/05 19:42:13 jim Exp $";

/*
 * revision history:
 *
 * $Log: dpmta_test4.c,v $
 * Revision 1.1  1997/09/05 19:42:13  jim
 * Original distribution.
 *
 * Revision 2.9  1997/05/07  21:28:07  chumphre
 * cubectr changed to cellctr and other minor changes related to the interface
 *
 * Revision 2.8  1997/05/07  19:31:52  chumphre
 * implement a uniform interface for rectangular and ||-piped cells
 *
 * Revision 2.7  1997/04/10  17:19:35  wrankin
 * added elapse timing code
 *
 * Revision 2.6  1997/03/12  20:15:17  wrankin
 * updates to timer codes
 *
 * Revision 2.5  1997/02/26  16:54:36  wrankin
 * added support for serial non-pvm compilation of DPMTA
 *
 * Revision 2.4  1996/11/18  19:29:47  wrankin
 * update to virial code to use multipole virial calculation
 * old virial code is turned on by defining OLDVIRIAL
 * direct test program now computes absolute virial tensor.
 *
 * Revision 2.3  1996/10/29  19:34:23  wrankin
 * changes to makefile to support other platforms
 * fix for multi-master communication interface
 * new global dumpdata routine to support direct pbc/macro verification
 *
 * Revision 2.2  1996/10/18  17:05:10  wrankin
 * improved PVM message passing structures.
 * made current and cleaned up T3D library code
 * added additional test codes
 * implement new heirarchical make structure
 *
 * Revision 2.1  1996/09/24  18:44:27  wrankin
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
 *
*/

/* include files */
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
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

main( int argc, char *argv[] ) {

   int i, j, k;                  /* loop counters */
   int info;                     /* return value info */
   int mytid;                    /* master's task id */
   int iter;                     /* number of iterations */
   int tids[64];                 /* tids of pmta slaves */
   int num_parts;                /* number of parts to simulate */
   PmtaInitData initdata;        /* DPMTA initialization data */
   PmtaParticlePtr particlelist; /* array of particle information */
   PmtaPartInfoPtr forcelist;    /* array of resulting force information */
   PmtaPartInfoPtr flist_lj;     /* array of resulting LJ force info */

#if defined VIRIAL || defined OLDVIRIAL
   double vpot;
   PmtaVector vfor;
   double vpot_lj;
   PmtaVector vfor_lj;
#endif

   /* check and load command line arguments */
   if (argc != 9) {
      fprintf(stderr,
         "%s <#procs> <#lvls> <pdb-file> <fft> <mp> <theta> <iter> <pbc>\n",
         argv[0]);
      exit(-1);
      }

   initdata.nprocs = atoi(argv[1]);
   initdata.nlevels = atoi(argv[2]);
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
   initdata.calling_num = 0;
   initdata.calling_tids = (int *)NULL;

   /* read in the contents of the PDB file */
   info = Read_PDB(argv[3], &num_parts, &(initdata.v1), &(initdata.v2),
                   &(initdata.v3), &(initdata.cellctr), &particlelist);
   if ( info < 0 ) {
      exit(-1);
   }

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

   forcelist = (PmtaPartInfoPtr)malloc(num_parts * sizeof(PmtaPartInfo));
#ifdef COMP_LJ
   flist_lj = (PmtaPartInfoPtr)malloc(num_parts * sizeof(PmtaPartInfo));
#else
   flist_lj = NULL;
#endif

#ifdef TIME
      gettimeofday(&trun,0);
      elapse_time = (double)trun.tv_sec +
	 ((double)trun.tv_usec/1000000);
#endif


   /*
   *  repeat for multiple iterations
   */

   for (i=0; i<iter; i++) { 

      /*
      *  send particles out to slaves and wait for results
      */

      PMTAforce(num_parts, particlelist, forcelist, flist_lj);

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


/****************************************************************
 *
 *  Read_PDB - read the contents of a PDB file containing all
 *    water and return the values in the pmrticle list array.
 *
 *  Taken from the PME code.
 */

int Read_PDB( char *filename,
	      int *nparts,
	      PmtaVector *v1,
	      PmtaVector *v2,
	      PmtaVector *v3,
	      PmtaVector *ctr,
	      PmtaParticlePtr *plist )
{

   int i,j,k;
   char s1[80];
   double cgh, cgo;
   FILE *infile;

   /* open the small.pdb file that contains box info and particle data */
   if ( (infile = fopen(filename,"r") ) == NULL) {
      fprintf(stderr,"Error: Unable to open PDB file %s\n",filename);
      return(-1); 
   }

   /* end open start reading (below) */
   fscanf(infile,"%s %s %s",s1,s1,s1); /* will read the text */
   fscanf(infile,"%lf",&(v1->x));
   fscanf(infile,"%d",nparts);

   /*
    * used for water system, hydrogen and oxygen charges
    * taken from Sigma
    */

   cgh = 7.471271;
   cgo = cgh * -2.0;

   /* we are assuming a cubic box here with the origin at (0,0,0)*/
   v2->y = v1->x;
   v3->z = v1->x;
   ctr->x = v1->x / 2.0;
   ctr->y = v2->y / 2.0;
   ctr->z = v3->z / 2.0;

   /* allocate main particle info array (x,y,z,cg)*/
   *plist = (PmtaParticlePtr)malloc( (*nparts) * sizeof(PmtaParticle) );
   if ( (*plist) == NULL ) {
      fprintf(stderr,"Error: could not alloc particle list\n");
      return(-1);
   }

   /* read from file8 small.pdb */
   for (i = 0; i < (*nparts); i+=3) {
      fscanf(infile, "%lf %lf %lf",
	     &(*plist)[i].p.x, &(*plist)[i].p.y, &(*plist)[i].p.z);
      fscanf(infile, "%lf %lf %lf",
	     &(*plist)[i+1].p.x, &(*plist)[i+1].p.y, &(*plist)[i+1].p.z);
      fscanf(infile, "%lf %lf %lf",
	     &(*plist)[i+2].p.x, &(*plist)[i+2].p.y, &(*plist)[i+2].p.z); 

      (*plist)[i].q = cgo;
      (*plist)[i+1].q = cgh;
      (*plist)[i+2].q = cgh;


#ifdef COMP_LJ
      /* eventually want to get the actual parameters here */
      (*plist)[i].a = 0.0;
      (*plist)[i].b = 0.0;
      (*plist)[i+1].a = 0.0;
      (*plist)[i+1].b = 0.0;
      (*plist)[i+2].a = 0.0;
      (*plist)[i+2].b = 0.0;
#endif

   } /* for i */

   fclose(infile);

} /* Read_PDB */
















