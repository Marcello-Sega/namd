/* 
*  dpmta_test3.c - PVM process for DPMTA implentation
*     using the distributed calling scheme and an SPMD
*     programming model.
*
*  w. t. rankin
*
*  Copyright (c) 1995 Duke University
*  All rights reserved
*
*
*  this program will instantiate and initialize the slave processes,
*  spawn multiple copies of an itself process, each of which 
*  will send particle data, and collect the resulting particle information
*  for the parallel FMA algorithm implementation.  The data will be returned
*  to the master (initial) process and dumped into a file.
*
*  this program is an example of how to use the PMTA interface routines
*  instantiating and calling the PMTA slave processes from a distributed
*  application.  it is also an example of how to perform SPMD processing
*  under PVM on both a workstation cluster and the Cray-T3D.
*
*/

static char rcsid[]="$Id: dpmta_test3.c,v 1.3 1997/09/29 23:58:46 jim Exp $";

/*
 * revision history:
 *
 * $Log: dpmta_test3.c,v $
 * Revision 1.3  1997/09/29 23:58:46  jim
 * Incorporated changes from version 2.6.1 of DPMTA.
 *   - fixes for bad handling of empty/invalid multipoles when
 *     using large processor sets.
 *   - moved functions that provide data mapping to processors.  master
 *     and slave routines now call the same function in dpmta_distmisc.c
 * Also, switched pvmc.h back to pvm3.h.
 *
 * Revision 1.2  1997/09/12 22:56:37  jim
 * Modifications to work with converse pvm.
 *
 * Revision 1.1  1997/09/05 19:42:13  jim
 * Original distribution.
 *
 * Revision 2.13  1997/05/12  20:45:14  chumphre
 * New interface bug fix
 *
 * Revision 2.12  1997/05/07  21:28:04  chumphre
 * cubectr changed to cellctr and other minor changes related to the interface
 *
 * Revision 2.11  1997/04/28  15:57:26  wrankin
 * cleaned up indentation in code
 *
 * Revision 2.10  1997/01/13  22:12:34  wrankin
 * general cleanup for certain m4 definitions
 *
 * Revision 2.9  1996/11/18  19:29:45  wrankin
 * update to virial code to use multipole virial calculation
 * old virial code is turned on by defining OLDVIRIAL
 * direct test program now computes absolute virial tensor.
 *
 * Revision 2.8  1996/11/01  02:26:49  wrankin
 * modifications to support cray t3d compilation
 * version update for 2.5 release
 *
 * Revision 2.7  1996/10/29  19:34:21  wrankin
 * changes to makefile to support other platforms
 * fix for multi-master communication interface
 * new global dumpdata routine to support direct pbc/macro verification
 *
 * Revision 2.6  1996/09/24  18:44:20  wrankin
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
 * Revision 2.5  1996/08/09  15:31:20  chumphre
 * Implemented noncubic (rectangular) boxes.
 *
 * Revision 2.4  1996/06/12  19:29:14  chumphre
 * Updated PMTAinit and PMTAforce interface.
 *
 * Revision 2.3  1995/11/29  22:29:42  wrankin
 * addition of periodic boundary code
 * beginning work on virial code.
 *
 * Revision 2.2  1995/07/26  01:52:59  wrankin
 * updated Makefiles
 * slight performance improvement for direct calculations
 * clean up test code example
 * fixed t3dlib for new multipole code
 *
 * Revision 2.1  1995/06/13  04:26:25  wrankin
 * No updates.
 * Sources placed under CVS.
 *
 * Revision 1.1  1995/04/24  04:29:22  wrankin
 * Initial revision
 *
 *
*/

/* include files */
#include <stdio.h>
#include <stdlib.h>
#include "pvm3.h"
#include "dpmta.h"

main( int argc, char *argv[] )
{
   int mytid;                /* process's task id */
   int ptid;                 /* parent's task id */
   int mypid;                /* my process number */
   int nprocs;               /* number of slave processes */
   int info;                 /* pvm return info */

   int i, j, k;              /* loop counters */

   int num_parts;            /* total number of particles */
   PmtaVector v1;            /* cell vector 1 */
   PmtaVector v2;            /* cell vector 2 */
   PmtaVector v3;            /* cell vector 3 */

   int iter;                 /* number of iterations */
   int tids[64];             /* tids of pmta slaves */
   int mtids[64];            /* tids of calling processes */

   PmtaVector cubecenter;    /* coordinates of cube center */

   PmtaInitData    initdata;     /* DPMTA Initialization Data */
   PmtaParticlePtr particlelist; /* array of particle information */
   PmtaPartInfoPtr forcelist;    /* array of resulting force information */
   PmtaPartInfoPtr flist_lj;     /* array of resulting LJ force info */

   PmtaParticlePtr mast_plist;    /* array of particle information */
   PmtaPartInfoPtr mast_flist;    /* array of resulting force information */
   PmtaPartInfoPtr mast_flist_lj; /* array of resulting LJ forces */

#if defined VIRIAL || defined OLDVIRIAL
   double vpot;
   PmtaVector vfor;
   double vpot_lj;
   PmtaVector vfor_lj;
#endif

   /*
   *  enroll and see if we are the initial process
   */

   mytid = pvm_mytid();
#ifdef CRAY
   mypid = pvm_get_PE(mytid);
   ptid = pvm_gettid("",0);
#else
   mypid = pvm_joingroup("GlobalGroup");
   ptid = pvm_gettid("GlobalGroup",0);
#endif


   /*
   *  if we have the first process, then we need to take care
   *  of spawning (if needed) and distribution the processing
   *  parameters
   */

   if ( mypid == 0 ) {

      /* check and load command line arguments */
      if (argc != 9) {
         fprintf(stderr,
            "%s <#procs> <#lvls> <#parts> <fft> <mp> <theta> <iter> <pbc>\n",
            argv[0]);
         pvm_exit();
         exit(-1);
      }

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

      /*
      *  spawn other procs, wait for all other procs to join the
      *  group and then set contents of tid array
      */

#ifdef CRAY
      for ( i=0; i<nprocs; i++ ) {
         mtids[i] = pvm_gettid("",i);
      }
#else
      mtids[0] = mytid;
      if ( nprocs > 1 ) {
         info = pvm_spawn("dpmta_test3", NULL, 0, "", nprocs-1, &(mtids[1]));
      }
#endif


      /* set other constants */

      initdata.mp_lj = 4;
      initdata.fftblock = 4;
      initdata.v1.x = v1.x = 32.0;
      initdata.v1.y = v1.y = 0.0;
      initdata.v1.z = v1.z = 0.0;
      initdata.v2.x = v2.x = 0.0;
      initdata.v2.y = v2.y = 32.0;
      initdata.v2.z = v2.z = 0.0;
      initdata.v3.x = v3.x = 0.0;
      initdata.v3.y = v3.y = 0.0;
      initdata.v3.z = v3.z = 32.0;
      initdata.cellctr.x = cubecenter.x = 0.0;
      initdata.cellctr.y = cubecenter.y = 0.0;
      initdata.cellctr.z = cubecenter.z = 0.0;

      initdata.calling_num = nprocs;
      initdata.calling_tids = mtids; 

      /* adjust number of parts by the numbr of processes */
      num_parts = num_parts / nprocs;

#ifdef DATADUMP
      /* allocate particle arrays */
      mast_plist =
         (PmtaParticlePtr)malloc(num_parts*nprocs*sizeof(PmtaParticle));
      if ( mast_plist == (PmtaParticlePtr)NULL ) {
         fprintf(stderr,"Error: Particle Malloc failed\n");
         exit(-1);
      }

      mast_flist =
         (PmtaPartInfoPtr)malloc(num_parts*nprocs*sizeof(PmtaPartInfo));
      if ( mast_flist == (PmtaPartInfoPtr)NULL ) {
         fprintf(stderr,"Error: Flist Malloc failed\n");
         exit(-1);
      }

#ifdef COMP_LJ
      mast_flist_lj =
         (PmtaPartInfoPtr)malloc(num_parts*nprocs*sizeof(PmtaPartInfo));
      if ( mast_flist_lj == (PmtaPartInfoPtr)NULL ) {
         fprintf(stderr,"Error: Flist_LJ Malloc failed\n");
         exit(-1);
      }
#endif
#endif

      /* send other processes the number of particles */
      for ( i=1; i<nprocs; i++ ) {
         pvm_initsend(PvmDataDefault);
         pvm_pkint(&iter,1,1);         
         pvm_pkint(&num_parts,1,1);
         pvm_pkdouble(&v1.x,3,1);
         pvm_pkdouble(&v2.x,3,1);
         pvm_pkdouble(&v3.x,3,1);
         pvm_pkdouble(&cubecenter.x,1,1);
         pvm_pkdouble(&cubecenter.y,1,1);
         pvm_pkdouble(&cubecenter.z,1,1);
         pvm_send(mtids[i],1);
      } /* for i */

      /* initialize slave tasks */
      PMTAinit(&(initdata), tids);

   } /* if mypid=0 */

   else {  /* not initial process */

      /* get initial message from master */
      pvm_recv(ptid,1);
      pvm_upkint(&iter,1,1);
      pvm_upkint(&num_parts,1,1);
      pvm_upkdouble(&v1.x,3,1);
      pvm_upkdouble(&v2.x,3,1);
      pvm_upkdouble(&v3.x,3,1);
      pvm_upkdouble(&cubecenter.x,1,1);
      pvm_upkdouble(&cubecenter.y,1,1);
      pvm_upkdouble(&cubecenter.z,1,1);

   } /* else */


   /* now all calling processes need to register */
   PMTAregister();


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
   if ( particlelist == (PmtaParticlePtr)NULL ) {
      fprintf(stderr,"Error: Particle Malloc failed\n");
      exit(-1);
   }

   forcelist = (PmtaPartInfoPtr)malloc(num_parts * sizeof(PmtaPartInfo));
   if ( particlelist == (PmtaParticlePtr)NULL ) {
      fprintf(stderr,"Error: Flist Malloc failed\n");
      exit(-1);
   }

   flist_lj = (PmtaPartInfoPtr)malloc(num_parts * sizeof(PmtaPartInfo));
   if ( particlelist == (PmtaParticlePtr)NULL ) {
      fprintf(stderr,"Error: Flist_LJ Malloc failed\n");
      exit(-1);
   }


   /* create all particles */

   srand48((long)mypid);

   for (i=0; i<num_parts; i++) {
      particlelist[i].p.x = (drand48() - 0.5) * v1.x + cubecenter.x;
      particlelist[i].p.y = (drand48() - 0.5) * v2.y + cubecenter.y;
      particlelist[i].p.z = (drand48() - 0.5) * v3.z + cubecenter.z;
      if ( (i%2) == 0 )
         particlelist[i].q = 1.0;
      else
         particlelist[i].q = -1.0;

      particlelist[i].a = 1.0e-3;
      particlelist[i].b = 1.0e-6;

      forcelist[i].f.x = 0.0;
      forcelist[i].f.y = 0.0;
      forcelist[i].f.z = 0.0;
      forcelist[i].v = 0.0;

      flist_lj[i].f.x = 0.0;
      flist_lj[i].f.y = 0.0;
      flist_lj[i].f.z = 0.0;
      flist_lj[i].v = 0.0;
   } /* for i */

   /*
   *  send particles out to slaves and wait for results
   *  repeat for multiple iterations
   */

   for (i=0; i<iter; i++) {
      PMTAforce(num_parts, particlelist, forcelist, flist_lj);

#if defined VIRIAL || defined OLDVIRIAL
      if ( mypid == 0 ) {
	 PMTAvirial( &vpot, &vfor, &vpot_lj, &vfor_lj );
	 fprintf(stderr,"Virial = %lf (%lf,%lf,%lf)\n",
		 vpot, vfor.x, vfor.y, vfor.z);
#ifdef COMP_LJ
	 fprintf(stderr,"LJ Virial = %lf (%lf,%lf,%lf)\n",
		 vpot_lj, vfor_lj.x, vfor_lj.y, vfor_lj.z);
#endif
      } /* if mypid */
#endif
   } /* for i */

   /*
   *  send response back to master
   */

#ifdef DATADUMP
   pvm_initsend(PvmDataInPlace);
   pvm_pkdouble(&(particlelist[0].p.x),num_parts*6,1);
   pvm_pkdouble(&(forcelist[0].f.x),num_parts*4,1);
#ifdef COMP_LJ
   pvm_pkdouble(&(flist_lj[0].f.x),num_parts*4,1);
#endif
   pvm_send(ptid,2);
#endif

   /*
   *  if we are the master, unpack particle data and print out
   */
   if ( mypid == 0 ) {
#ifdef DATADUMP
      for ( i=0; i<nprocs; i++ ) {
         pvm_recv(-1,2);
         pvm_upkdouble(&(mast_plist[i*num_parts].p.x),
            num_parts*6,1);
         pvm_upkdouble(&(mast_flist[i*num_parts].f.x),
            num_parts*4,1);
#ifdef COMP_LJ
         pvm_upkdouble(&(mast_flist_lj[i*num_parts].f.x),
            num_parts*4,1);
#endif
      } /* for i */

      Dump_Results((num_parts*nprocs), &(initdata), mast_plist,
		   mast_flist, mast_flist_lj);

#endif

      PMTAexit();


   } /* if mypid */


   /*
   *  clean up and kill off slaves
   */

   pvm_exit();
   exit(0);

} /* main */

