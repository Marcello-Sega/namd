/* 
*  dpmta_direct.c - PVM process that computes the electrostatics
*     using the naive N^2 direct method but using an SPMD
*     programming model.
*
*  w. t. rankin
*
*  Copyright (c) 1995 Duke University
*  All rights reserved
*
*
*  this program will read in a file containing particle position data,
*  instantiate and initialize the slave processes,
*  spawn multiple copies of itself, each of which will receive the
*  entire particle data set, and compute the portion corrresponding
*  to the processor id.  the data will be returned
*  to the master (initial) process and dumped into a file.
*
*  this program is intended to be used to verify the output of the
*  various dpmta test programs.
*
*
*/

static char rcsid[]="$Id: dpmta_direct.c,v 1.2 1997/09/12 22:56:29 jim Exp $";

/*
 * revision history:
 *
 * $Log: dpmta_direct.c,v $
 * Revision 1.2  1997/09/12 22:56:29  jim
 * Modifications to work with converse pvm.
 *
 * Revision 1.1  1997/09/05 19:41:51  jim
 * Original distribution.
 *
 * Revision 2.10  1997/05/07  18:59:16  chumphre
 * major changes to implement ||-piped cells
 *
 * Revision 2.9  1997/03/26  17:45:59  wrankin
 * fixed mis-type to allow serial execution
 *
 * Revision 2.8  1997/03/07  00:46:12  wrankin
 * converted dpmta_direct to work serial without pvm
 *
 * Revision 2.7  1996/11/18  19:29:24  wrankin
 * update to virial code to use multipole virial calculation
 * old virial code is turned on by defining OLDVIRIAL
 * direct test program now computes absolute virial tensor.
 *
 * Revision 2.6  1996/10/28  23:01:52  wrankin
 * additions to test routines to provide processing paramenters as
 *   part of output file.
 * dpmta_direct will now perform macroscopic computations, reading
 *   in processing parameters from particle position file.
 *
 * Revision 2.5  1996/09/24  18:41:40  wrankin
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
 * Revision 2.4  1996/08/14  16:06:58  wrankin
 * Fixed LJ computations
 *
 * Revision 2.3  1996/02/29  21:13:21  wrankin
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
 * Revision 2.2  1995/10/01  21:45:45  wrankin
 * implemented LJ multipole calculations - compile time option only
 *
 * Revision 2.1  1995/06/13  04:26:01  wrankin
 * No updates.
 * Sources placed under CVS.
 *
 * Revision 1.1  1995/04/24  04:17:35  wrankin
 * Initial revision
 *
 *
*/

/* include files */
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#ifndef SERIAL
#include "pvmc.h"
#endif
#include "dpmta.h"

/* prototype me baby */
int MAC( double, double, double, double );

main( int argc, char *argv[] )
{
   int mytid;                /* process's task id */
   int ptid;                 /* parent's task id */
   int mypid;                /* my process number */
   int nprocs;               /* number of slave processes */
   int nparts;               /* number of particles */

   int i, j, k;              /* loop counters */
   int s_indx, e_indx;       /* loop index limits */

   int tids[64];             /* tids of pmta slaves */

   PmtaParticlePtr plist;    /* array of particle information */
   PmtaPartInfoPtr flist;    /* array of resulting force information */

   int pbc, kterm;           /* multipole pbc and macro flags */
   double theta;             /* multipole sep criteria */
#ifndef PIPED
   PmtaVector cl, cc;        /* cell edge lengths and center */
#else
   PmtaVector v1, v2, v3, cc;/* ||-Piped vectors and center */
#endif
   char ctmp[80];            /* temporary string for reading in data */

   int x, y, z;              /* pbc loop counters */
   int xmax, ymax, zmax;     /* pbc boundaries */
   double rad1, rad2;        /* radius used in MAC */
   PmtaVector sep;           /* cell edge lengths and center */

   double wtq;
   double dx, dy, dz;        /* processing parameters */
   double ir, ir2;           /* processing parameters */
   double ir_q, ir3_q;
   double ir3dx, ir3dy, ir3dz;

   FILE *fpin, *fpout;       /* input and output file pointers */

#ifdef COMP_LJ
   double wta, wtb;
   double ir6, ir6_a, ir8_a, ir12_b, ir14_b, ir_lj;
   PmtaPartInfoPtr flist_lj;
#endif

#if defined VIRIAL || defined OLDVIRIAL
   PmtaVector virf;          /* virial tensor */
   double virp;              /* virial potential */
   double ftmp[4];
#endif

#ifdef MACROSCOPIC
   int x2, y2, z2;           /* macro loop counters */
   int x3, y3, z3;           /* macro loop counters */
   int klvl;                 /* macroscopic level */
   int ncell, kdist;         /* some distance counters */
   PmtaVector *mv;           /* list of macroscopic interactions */
   PmtaVector sep2;          /* cell edge lengths and center */
#endif


#ifdef SERIAL
   mypid = 0;
#else
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

#endif

   /*
   *  if we have the first process, then we need to take care
   *  of spawning (if needed) reading in the input data and 
   *  distribution the processing
   *  parameters
   */

   if ( mypid == 0 ) {

#ifdef COMP_LJ
      /* check and load command line arguments */
      if (argc != 5) {
         fprintf(stderr,
            "%s <#procs> <infile> <C-outfile> <LJ-outfile>\n",
            argv[0]);
#ifndef SERIAL
         pvm_exit();
#endif
         exit(-1);
         }
#else
      /* check and load command line arguments */
      if (argc != 4) {
         fprintf(stderr,
            "%s <#procs> <infile> <outfile>\n",
            argv[0]);
#ifndef SERIAL
         pvm_exit();
#endif
         exit(-1);
         }
#endif

#ifdef SERIAL
      nprocs = 1;
#else
      nprocs = atoi(argv[1]);
#endif

      /*
      *  spawn other procs, wait for all other procs to join the
      *  group and then set contents of tid array
      */

#ifndef SERIAL
#ifdef CRAY
      for ( i=0; i<nprocs; i++ ) {
         tids[i] = pvm_gettid("",i);
      }
#else
      tids[0] = mytid;
      if ( nprocs > 1 ) {
	 pvm_spawn("dpmta_direct", NULL, 0, "", nprocs-1, &(tids[1]));
      }
#endif
#endif

      /* read in the particle data */
      fpin = fopen(argv[2],"r");

#ifndef PIPED
      fscanf(fpin, "%s %lg %lg %lg\n", ctmp, &(cl.x), &(cl.y), &(cl.z));
#else
      fscanf(fpin, "%s %lg %lg %lg\n", ctmp, &(v1.x), &(v1.y), &(v1.z));
      fscanf(fpin, "%s %lg %lg %lg\n", ctmp, &(v2.x), &(v2.y), &(v2.z));
      fscanf(fpin, "%s %lg %lg %lg\n", ctmp, &(v3.x), &(v3.y), &(v3.z));
#endif
      fscanf(fpin, "%s %lg %lg %lg\n", ctmp, &(cc.x), &(cc.y), &(cc.z));
      fscanf(fpin, "%s %d %d\n", ctmp, &pbc, &kterm);
      fscanf(fpin, "%s %lg\n", ctmp, &theta);
      fscanf(fpin, "%s %d\n", ctmp, &nparts);

#ifndef PIPED
      fprintf(stderr, "cell length:  %lg %lg %lg\n", cl.x, cl.y, cl.z);
#else
      fprintf(stderr, "||-Piped vector 1: %lg %lg %lg\n", v1.x, v1.y, v1.z);
      fprintf(stderr, "||-Piped vector 2: %lg %lg %lg\n", v2.x, v2.y, v2.z);
      fprintf(stderr, "||-Piped vector 3: %lg %lg %lg\n", v3.x, v3.y, v3.z);
#endif
      fprintf(stderr, "cell center: %lg %lg %lg\n", cc.x, cc.y, cc.z);
      fprintf(stderr, "pbc/kterm: %d %d\n", pbc, kterm);
      fprintf(stderr, "theta: %lg\n", theta);
      fprintf(stderr, "nparts: %d\n", nparts);

      /* allocate particle arrays */
      plist = (PmtaParticlePtr)malloc(nparts*sizeof(PmtaParticle));

      for ( i=0; i<nparts; i++ ) {
#ifdef COMP_LJ
         fscanf(fpin, "%lf %lf %lf %lf %lf %lf", 
                &(plist[i].p.x), &(plist[i].p.y), &(plist[i].p.z),
		&(plist[i].q), &(plist[i].a), &(plist[i].b));
#else
         fscanf(fpin, "%lf %lf %lf %lf ",
		&(plist[i].p.x), &(plist[i].p.y), &(plist[i].p.z),
		&(plist[i].q)); 
#endif
      }
      fclose(fpin);

#ifndef SERIAL
      /* send other processes the particle data */
      for ( i=1; i<nprocs; i++ ) {
         pvm_initsend(PvmDataRaw);
         pvm_pkint(&nprocs,1,1);
#ifndef PIPED
         pvm_pkdouble(&(cl.x),3,1);
#else
         pvm_pkdouble(&(v1.x),3,1);
         pvm_pkdouble(&(v2.x),3,1);
         pvm_pkdouble(&(v3.x),3,1);
#endif
         pvm_pkdouble(&(cc.x),3,1);
         pvm_pkdouble(&theta,1,1);
         pvm_pkint(&pbc,1,1);
         pvm_pkint(&kterm,1,1);
         pvm_pkint(&nparts,1,1);
         pvm_pkdouble(&(plist[0].p.x),nparts*6,1);
         pvm_send(tids[i],1);
      }
#endif

   } /* if mypid */

   /****************************************************************/

#ifndef SERIAL
   else {  /* not initial process */

      /* get initial message from master */
      pvm_recv(ptid,1);
      pvm_upkint(&nprocs,1,1);
#ifndef PIPED
      pvm_upkdouble(&(cl.x),3,1);
#else
      pvm_upkdouble(&(v1.x),3,1);
      pvm_upkdouble(&(v2.x),3,1);
      pvm_upkdouble(&(v3.x),3,1);
#endif
      pvm_upkdouble(&(cc.x),3,1);
      pvm_upkdouble(&theta,1,1);
      pvm_upkint(&pbc,1,1);
      pvm_upkint(&kterm,1,1);

      pvm_upkint(&nparts,1,1);

      /* allocate particle arrays and unpack data */
      plist = (PmtaParticlePtr)malloc(nparts*sizeof(PmtaParticle));
      pvm_upkdouble(&(plist[0].p.x),nparts*6,1);
   }
#endif

   /*
   *  at this point, all processors now have all particle data.
   *  we need to initialze the force lists and them compute which
   *  sections of data that we will need to work on.
   */

   flist = (PmtaPartInfoPtr)malloc(nparts*sizeof(PmtaPartInfo));
   for ( i=0; i<nparts; i++ ) {
      flist[i].f.x = 0.0;
      flist[i].f.y = 0.0;
      flist[i].f.z = 0.0;
      flist[i].v = 0.0;
   }

#ifdef COMP_LJ
   flist_lj = (PmtaPartInfoPtr)malloc(nparts*sizeof(PmtaPartInfo));
   for ( i=0; i<nparts; i++ ) {
      flist_lj[i].f.x = 0.0;
      flist_lj[i].f.y = 0.0;
      flist_lj[i].f.z = 0.0;
      flist_lj[i].v = 0.0;
   }
#endif

#if defined VIRIAL || defined OLDVIRIAL
   virf.x = 0.0;
   virf.y = 0.0;
   virf.z = 0.0;
   virp = 0.0;
#endif

   s_indx = (nparts/nprocs)*mypid;
   if ( mypid == (nprocs-1) )
      e_indx = nparts;
   else
      e_indx = (nparts/nprocs)*(mypid+1);


   for ( i=s_indx; i<e_indx; i++ ) {
      for ( j=0; j<nparts; j++ ) {
         if ( i != j ) {

            wtq = plist[i].q * plist[j].q;

            dx = plist[j].p.x - plist[i].p.x;
            dy = plist[j].p.y - plist[i].p.y;
            dz = plist[j].p.z - plist[i].p.z;

            ir2 = 1.0/(dx*dx + dy*dy + dz*dz);
	    ir = sqrt(ir2);

            ir_q = wtq * ir;
            flist[i].v += ir_q;

            ir3_q = ir_q * ir2;

            ir3dx = ir3_q * dx;
	    ir3dy = ir3_q * dy;
	    ir3dz = ir3_q * dz;
	    flist[i].f.x -= ir3dx;
	    flist[i].f.y -= ir3dy;
	    flist[i].f.z -= ir3dz;

#if defined VIRIAL || defined OLDVIRIAL
            virp -= ir_q * 0.5;
            virf.x -= ir3dx * dx * 0.5;
            virf.y -= ir3dy * dy * 0.5;
            virf.z -= ir3dz * dz * 0.5;
#endif

#ifdef COMP_LJ
            wta = plist[i].a * plist[j].a;
            wtb = plist[i].b * plist[j].b;

            ir6 = ir2 * ir2 *ir2;
            ir6_a = ir6 * wta;
     	    ir8_a = ir6_a * ir2;
	    ir12_b = ir6 * ir6 * wtb;
	    ir14_b = ir12_b * ir2;

	    flist_lj[i].v += ir12_b - ir6_a;
	    ir_lj = 6.0 * ( 2.0 * ir14_b - ir8_a );
	    flist_lj[i].f.x -= ir_lj * dx;
	    flist_lj[i].f.y -= ir_lj * dy;
	    flist_lj[i].f.z -= ir_lj * dz;

#endif

         } /* j != i */
      } /* for j */
   } /* for i */

   /*
    *  compute pbc's if needed 
    */

   if ( pbc == 1 ) {
#ifdef PIPED
      fprintf(stderr, "PBC's not implemented for ||-pipeds.\n");
      fprintf(stderr, "Results will not include PBC calculations.\n");
#else
      rad1 = sqrt(cl.x*cl.x + cl.y*cl.y + cl.z*cl.z) / 2.0;
      xmax = 0;
      while ( MAC(rad1, rad1, (double)(xmax)*cl.x, theta) == 0 )
	 xmax++;
      ymax = 0;
      while ( MAC(rad1, rad1, (double)(ymax)*cl.y, theta) == 0 )
	 ymax++;
      zmax = 0;
      while ( MAC(rad1, rad1, (double)(zmax)*cl.z, theta) == 0 )
	 zmax++;

      for ( x=-xmax; x<=xmax; x++ ) {
	 sep.x = (double)x * cl.x;
	 for ( y=-ymax; y<=ymax; y++ ) {
	    sep.y = (double)y * cl.y;
	    for ( z=-zmax; z<=zmax; z++ ) {
 	       sep.z = (double)z * cl.z;
	       if (( x != 0 )||( y != 0 )||( z != 0 )) {
		  rad2 = sqrt(sep.x*sep.x + sep.y*sep.y + sep.z*sep.z);
		  if ( MAC(rad1, rad1, rad2, theta) == 0 ) {
		     for ( i=s_indx; i<e_indx; i++ ) {
		        for ( j=0; j<nparts; j++ ) {

			   wtq = plist[i].q * plist[j].q;

			   dx = plist[j].p.x - plist[i].p.x + sep.x;
			   dy = plist[j].p.y - plist[i].p.y + sep.y;
			   dz = plist[j].p.z - plist[i].p.z + sep.z;

			   ir2 = 1.0/(dx*dx + dy*dy + dz*dz);
			   ir = sqrt(ir2);

			   ir_q = wtq * ir;
			   flist[i].v += ir_q;

			   ir3_q = ir_q * ir2;

			   ir3dx = ir3_q * dx;
			   ir3dy = ir3_q * dy;
			   ir3dz = ir3_q * dz;
			   flist[i].f.x -= ir3dx;
			   flist[i].f.y -= ir3dy;
			   flist[i].f.z -= ir3dz;

#if defined VIRIAL || defined OLDVIRIAL
			   virp -= ir_q * 0.5;
			   virf.x -= ir3dx * dx * 0.5;
			   virf.y -= ir3dy * dy * 0.5;
			   virf.z -= ir3dz * dz * 0.5;
#endif

#ifdef COMP_LJ
			   wta = plist[i].a * plist[j].a;
			   wtb = plist[i].b * plist[j].b;

			   ir6 = ir2 * ir2 *ir2;
			   ir6_a = ir6 * wta;
			   ir8_a = ir6_a * ir2;
			   ir12_b = ir6 * ir6 * wtb;
			   ir14_b = ir12_b * ir2;

			   flist_lj[i].v += ir12_b - ir6_a;
			   ir_lj = 6.0 * ( 2.0 * ir14_b - ir8_a );
			   flist_lj[i].f.x -= ir_lj * dx;
			   flist_lj[i].f.y -= ir_lj * dy;
			   flist_lj[i].f.z -= ir_lj * dz;
#endif
			} /* for j */
		     } /* for i */
		  } /* if MAC */
	       } /* if not 0 */
	    } /* for z */
	 } /* for y */
      } /* for z */

#endif /* PIPED */

   } /* if pbc */

#ifdef MACROSCOPIC

   /*
    *  compute macroscopic expansions if needed 
    */

   if ( kterm > 0 ) {

      for ( x=-xmax; x<=xmax; x++ ) {
	 sep.x = (double)x * cl.x;
	 for ( y=-ymax; y<=ymax; y++ ) {
	    sep.y = (double)y * cl.y;
	    for ( z=-zmax; z<=zmax; z++ ) {
 	       sep.z = (double)z * cl.z;
	       rad2 = sqrt(sep.x*sep.x + sep.y*sep.y + sep.z*sep.z);
	       if ( MAC(rad1, rad1, rad2, theta) == 0 ) {
		  for ( x2=-1; x2<=1; x2++ ) {
		     sep2.x = (double)(3*x+x2) * cl.x;
		     for ( y2=-1; y2<=1; y2++ ) {
			sep2.y = (double)(3*y+y2) * cl.y;
			for ( z2=-1; z2<=1; z2++ ) {
			   sep2.z = (double)(3*z+z2) * cl.z;
			   rad2 = sqrt(sep2.x*sep2.x + sep2.y*sep2.y
				       + sep2.z*sep2.z);
			   if ( MAC(rad1, rad1, rad2, theta) == 1 ) {
			      kdist = 1;
			      for ( klvl=0; klvl<kterm; klvl++ ) {
				 ncell = kdist/2;
				 for ( x3=(-ncell); x3<=ncell; x3++ ) {
				    for ( y3=(-ncell); y3<=ncell; y3++ ) {
				       for ( z3=(-ncell); z3<=ncell; z3++ ) {

					  for ( i=s_indx; i<e_indx; i++ ) {
					     for ( j=0; j<nparts; j++ ) {

						wtq = plist[i].q * plist[j].q;

						dx = plist[j].p.x - plist[i].p.x + (double)(kdist*(3*x+x2)+x3)*cl.x;
						dy = plist[j].p.y - plist[i].p.y + (double)(kdist*(3*y+y2)+y3)*cl.y;
						dz = plist[j].p.z - plist[i].p.z + (double)(kdist*(3*z+z2)+z3)*cl.z;

						ir2 = 1.0/(dx*dx + dy*dy + dz*dz);
						ir = sqrt(ir2);

						ir_q = wtq * ir;
						flist[i].v += ir_q;

						ir3_q = ir_q * ir2;

						ir3dx = ir3_q * dx;
						ir3dy = ir3_q * dy;
						ir3dz = ir3_q * dz;
						flist[i].f.x -= ir3dx;
						flist[i].f.y -= ir3dy;
						flist[i].f.z -= ir3dz;

#if defined VIRIAL || defined OLDVIRIAL
						virp -= ir_q * 0.5;
						virf.x -= ir3dx * dx * 0.5;
						virf.y -= ir3dy * dy * 0.5;
						virf.z -= ir3dz * dz * 0.5;
#endif

#ifdef COMP_LJ
						wta = plist[i].a * plist[j].a;
						wtb = plist[i].b * plist[j].b;

						ir6 = ir2 * ir2 *ir2;
						ir6_a = ir6 * wta;
						ir8_a = ir6_a * ir2;
						ir12_b = ir6 * ir6 * wtb;
						ir14_b = ir12_b * ir2;

						flist_lj[i].v += ir12_b - ir6_a;

						ir_lj = 6.0 * ( 2.0 * ir14_b - ir8_a );
						flist_lj[i].f.x -= ir_lj * dx;
						flist_lj[i].f.y -= ir_lj * dy;
						flist_lj[i].f.z -= ir_lj * dz;
#endif
					     } /* for j */
					  } /* for i */
				       } /* for z3 */
				    } /* for y3 */
				 } /* for x3 */
				 kdist *= 3;
			      } /* for klvl */
			   } /* if MAC */
			} /* for z2 */
		     } /* for y2 */
		  } /* for x2 */
	       } /* if MAC */
	    } /* for z */
	 } /* for y */
      } /* for z */

   } /* if kterm */

#endif

   /*
   *  if we are the master process, collect other results
   *  otherwise, return particles to the first process
   */

   if ( mypid == 0 ) {

#ifndef SERIAL
      /*
      *  collect results from slaves
      */

      for ( i=1; i<nprocs; i++ ) {
         pvm_recv(-1,2);
         pvm_upkint(&(j),1,1);
         pvm_upkint(&(k),1,1);
         pvm_upkdouble(&(flist[k].f.x),j*4,1);
#ifdef COMP_LJ
         pvm_upkdouble(&(flist_lj[k].f.x),j*4,1);
#endif
#if defined VIRIAL || defined OLDVIRIAL
         pvm_upkdouble(ftmp,4,1);
         virp += ftmp[0];
         virf.x += ftmp[1];
         virf.y += ftmp[2];
         virf.z += ftmp[3];
#endif
      }
      
#endif

      /*
      *  print results
      */
      fpout = fopen(argv[3],"w");
      for ( i=0; i<nparts; i++ ) {
         fprintf(fpout,"%20.16lf %20.16lf %20.16lf %20.16lf\n",
            flist[i].f.x, flist[i].f.y, flist[i].f.z, flist[i].v);
      }
      fclose(fpout);

#ifdef COMP_LJ
      fpout = fopen(argv[4],"w");
      for ( i=0; i<nparts; i++ ) {
         fprintf(fpout,"%20.16lf %20.16lf %20.16lf %20.16lf\n",
            flist_lj[i].f.x, flist_lj[i].f.y, flist_lj[i].f.z, flist_lj[i].v);
      }
      fclose(fpout);
#endif

#if defined VIRIAL || defined OLDVIRIAL
      fprintf(stderr,"Virial = %lf (%lf,%lf,%lf)\n",
	      virp, virf.x, virf.y, virf.z);
#endif

   } /* if mypid == 0 */


#ifndef SERIAL
   else {

      /*
      *  send response back to master
      */

      nparts = e_indx - s_indx;
      pvm_initsend(PvmDataRaw);

      pvm_pkint(&(nparts),1,1);
      pvm_pkint(&(s_indx),1,1);
      pvm_pkdouble(&(flist[s_indx].f.x),nparts*4,1);
#ifdef COMP_LJ
      pvm_pkdouble(&(flist_lj[s_indx].f.x),nparts*4,1);
#endif
#if defined VIRIAL || defined OLDVIRIAL
      pvm_pkdouble(&(virp),1,1);
      pvm_pkdouble(&(virf.x),3,1);
#endif
      pvm_send(ptid,2);
   } /* mypid != 0 */


   /*
   *  exit from pvm
   */
   pvm_exit();

#endif

   exit(0);
}


int MAC( double r1, double r2, double rsep, double theta )
{
   if ( (r1+r2) <= (theta*rsep) )
      return(1);
   else
      return(0);
}
