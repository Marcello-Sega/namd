/* 
*  dpmta_test2.c - master PVM process for DPMTA implentation
*     useing the distributed calling scheme.
*
*  w. t. rankin
*
*  Copyright (c) 1994 Duke University
*  All rights reserved
*
*  this program will instantiate and initialize the slave processes,
*  spawn multiple copies of an application slave process, each of which 
*  will send particle data, and collect the resulting particle information
*  for the parallel FMA algorithm implementation.  The data will be returned
*  to teh master process and dumped into a file.
*
*  this program is an example of how to use the PMTA interface routines
*  instantiating and calling the PMTA slave processes from a distributed
*  application.
*
*/

static char rcsid[]="$Id: dpmta_test2.c,v 1.2 1997/09/12 22:56:36 jim Exp $";

/*
 * revision history:
 *
 * $Log: dpmta_test2.c,v $
 * Revision 1.2  1997/09/12 22:56:36  jim
 * Modifications to work with converse pvm.
 *
 * Revision 1.1  1997/09/05 19:42:11  jim
 * Original distribution.
 *
 * Revision 2.6  1997/05/12  20:43:00  chumphre
 * New interface bug fix.
 *
 * Revision 2.5  1997/05/07  21:28:00  chumphre
 * cubectr changed to cellctr and other minor changes related to the interface
 *
 * Revision 2.4  1996/10/29  19:34:19  wrankin
 * changes to makefile to support other platforms
 * fix for multi-master communication interface
 * new global dumpdata routine to support direct pbc/macro verification
 *
 * Revision 2.3  1996/09/24  18:44:08  wrankin
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
 * Revision 2.2  1996/02/06  00:55:05  wrankin
 * update for version 2.3 release.
 *
 * Revision 2.1  1995/06/13  04:26:23  wrankin
 * No updates.
 * Sources placed under CVS.
 *
 * Revision 1.2  1995/02/24  21:23:32  wrankin
 * added interaction list parameter theta as a command line argument
 *
 * Revision 1.1  1994/11/30  17:50:21  wrankin
 * Initial revision
 *
 *
 *
*/

/* include files */
#include <stdio.h>
#include <stdlib.h>
#include "pvmc.h"
#include "dpmta.h"

/* defines */
#define MAXPROCS 64

main( int argc, char *argv[] )
{
   int mytid;                /* master's task id */
   int nprocs;               /* number of slave processes */

   int i, j, k;              /* loop counters */

   int num_parts;            /* total number of particles */
   int iter;                 /* number of iterations */

   int tids[MAXPROCS];       /* tids of pmta slaves */
   int mtids[MAXPROCS];      /* tids of pmta master copies */

   PmtaInitData initdata;        /* dpmta initialization data */
   PmtaParticlePtr particlelist; /* array of particle information */
   PmtaPartInfoPtr forcelist;    /* array of resulting force information */
   PmtaPartInfoPtr flist_lj;     /* array of resulting force information */


   /* check and load command line arguments */
   if (argc != 9) {
      fprintf(stderr,
         "%s <#procs> <#lvls> <#parts> <fft> <mp> <theta> <iter> <pbc>\n",
         argv[0]);
      exit(-1);
      }

   /* enroll master in pvm */
   mytid = pvm_mytid();


   initdata.nprocs = (nprocs = atoi(argv[1]));
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
   initdata.v1.x = 32.0;
   initdata.v1.y = 0.0;
   initdata.v1.z = 0.0;
   initdata.v2.x = 0.0;
   initdata.v2.y = 32.0;
   initdata.v2.z = 0.0;
   initdata.v3.x = 0.0;
   initdata.v3.y = 0.0;
   initdata.v3.z = 32.0;
   initdata.cellctr.x = 16.0;
   initdata.cellctr.y = 16.0;
   initdata.cellctr.z = 16.0;

   /* spawn other (sub)masters */
   pvm_spawn("dpmta_test2b", NULL, 0, "", nprocs, mtids);

   /* divide total number of particles by the number of procs */
   num_parts /= nprocs;

   /* send initialization message */
   for ( i=0; i<nprocs; i++ ) {
      pvm_initsend(PvmDataDefault);
      pvm_pkint(&iter,1,1);
      pvm_pkint(&num_parts,1,1);
      pvm_pkdouble(&(initdata.v1.x),1,1);
      pvm_pkdouble(&(initdata.v2.y),1,1);
      pvm_pkdouble(&(initdata.v3.z),1,1);
      pvm_pkdouble(&(initdata.cellctr.x),3,1);
      pvm_send(mtids[i],1);
      }

   initdata.calling_num = nprocs;
   initdata.calling_tids = mtids;

   /* initialize slave tasks */
   PMTAinit(&(initdata), tids);

   /* allocate particle arrays */
   particlelist =
      (PmtaParticlePtr)malloc(nprocs*num_parts*sizeof(PmtaParticle));
   forcelist =
      (PmtaPartInfoPtr)malloc(nprocs*num_parts*sizeof(PmtaPartInfo));
   flist_lj =
      (PmtaPartInfoPtr)malloc(nprocs*num_parts*sizeof(PmtaPartInfo));


   /* wait for return messages from slaves */
   for ( i=0; i<nprocs; i++ ) {
      pvm_recv(-1,2);
      pvm_upkdouble(&(particlelist[i*num_parts].p.x),
         num_parts*6,1);
      pvm_upkdouble(&(forcelist[i*num_parts].f.x),
         num_parts*4,1);
      pvm_upkdouble(&(flist_lj[i*num_parts].f.x),
         num_parts*4,1);
      }

#ifdef DATADUMP
   Dump_Results(nprocs*num_parts, &(initdata), particlelist, forcelist, flist_lj);
#endif

   /*
   *  clean up and kill off slaves
   */

   PMTAexit();

   pvm_exit();
   exit(0);
}

