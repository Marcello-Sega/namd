/* 
*  dpmta_test5.c - PVM process for DPMTA implentation
*     using the distributed calling scheme and an SPMD
*     programming model.
*
*  w. t. rankin
*
*  Copyright (c) 1997 Duke University
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

static char rcsid[]="$Id: dpmta_test5.c,v 1.2 1997/09/12 22:56:38 jim Exp $";

/*
 * revision history:
 *
 * $Log: dpmta_test5.c,v $
 * Revision 1.2  1997/09/12 22:56:38  jim
 * Modifications to work with converse pvm.
 *
 * Revision 1.1  1997/09/05 19:42:14  jim
 * Original distribution.
 *
 * Revision 2.2  1997/05/07  21:28:09  chumphre
 * cubectr changed to cellctr and other minor changes related to the interface
 *
 * Revision 2.1  1997/01/13  22:09:59  wrankin
 * added supprt for dpmta_test5.c which does MIMD with PDB file input
 *
 *
*/

/* include files */
#include <stdio.h>
#include <stdlib.h>
#include "pvmc.h"
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
   PmtaVector v1;            /* cell vector 1*/
   PmtaVector v2;            /* cell vector 2*/
   PmtaVector v3;            /* cell vector 3*/

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


      /* set other constants */

      initdata.mp_lj = 4;
      initdata.fftblock = 4;
      initdata.calling_num = nprocs;
      initdata.calling_tids = mtids; 

      /* read in the contents of the PDB file */
      info = Read_PDB(argv[3],&num_parts,&(v1),&(v2),&(v3),&(cubecenter),
		   &mast_plist);
      if ( info < 0 ) {
         exit(-1);
      }

      initdata.v1.x = v1.x;
      initdata.v2.y = v2.y;
      initdata.v3.z = v3.z;
      initdata.cellctr.x = cubecenter.x;
      initdata.cellctr.y = cubecenter.y;
      initdata.cellctr.z = cubecenter.z;


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
         info = pvm_spawn("dpmta_test5", NULL, 0, "", nprocs-1, &(mtids[1]));
      }
#endif

      /* allocate particle arrays */

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


      /* send other processes the number of particles */
      for ( i=1; i<nprocs; i++ ) {
         pvm_initsend(PvmDataDefault);
         pvm_pkint(&iter,1,1);         
         pvm_send(mtids[i],1);
      }

      /* initialize slave tasks */
      PMTAinit(&(initdata), tids);

   }

   else {  /* not initial process */


      /* get initial message from master */
      pvm_recv(ptid,1);
      pvm_upkint(&iter,1,1);

      num_parts = 0;
      mast_plist = (PmtaParticlePtr)NULL;
      mast_flist = (PmtaPartInfoPtr)NULL;
      mast_flist_lj = (PmtaPartInfoPtr)NULL;
   }


   /* now all calling processes need to register */
   PMTAregister();


   /*
   *  send particles out to slaves and wait for results
   *  repeat for multiple iterations
   */

   for (i=0; i<iter; i++) {
      PMTAforce(num_parts, mast_plist, mast_flist, mast_flist_lj);

#if defined VIRIAL || defined OLDVIRIAL
      if ( mypid == 0 ) {
	 PMTAvirial( &vpot, &vfor, &vpot_lj, &vfor_lj );
	 fprintf(stderr,"Virial = %lf (%lf,%lf,%lf)\n",
		 vpot, vfor.x, vfor.y, vfor.z);
#ifdef COMP_LJ
	 fprintf(stderr,"LJ Virial = %lf (%lf,%lf,%lf)\n",
		 vpot_lj, vfor_lj.x, vfor_lj.y, vfor_lj.z);
#endif
      }
#endif
   }

   /*
   *  send response back to master
   */

   pvm_initsend(PvmDataInPlace);
   pvm_send(ptid,2);

   /*
   *  if we are the master, unpack particle data and print out
   */
   if ( mypid == 0 ) {
      for ( i=0; i<nprocs; i++ ) {
         pvm_recv(-1,2);
      } /* for i */


#ifdef DATADUMP
      Dump_Results(num_parts, &(initdata), mast_plist,
		   mast_flist, mast_flist_lj);
#endif

      PMTAexit();


   } /* if mypid */

   /*
   *  clean up and kill off slaves
   */

   pvm_exit();
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
	      PmtaVector *box,
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
   fscanf(infile,"%lf",&(box->x));
   fscanf(infile,"%d",nparts);

   /* used for water system, hydrogen and oxygen charges */
   /* cgh = sqrt(332.17752) * .417; */
   /* cgo = cgh * -2.; */

   cgh = 0.410;
   cgo = -0.820;

   /* we are assuming a cubic box here with the origin at (0,0,0)*/
   box->y = box->x;
   box->z = box->x;
   ctr->x = box->x / 2.0;
   ctr->y = box->y / 2.0;
   ctr->z = box->z / 2.0;

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
