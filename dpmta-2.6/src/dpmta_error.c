/* 
*  dpmta_error.c - compute the error between multipole and
*    direct results.
*
*  w. t. rankin
*
*  Copyright (c) 1996 Duke University
*  All rights reserved
*
*/

static char rcsid[]="$Id: dpmta_error.c,v 1.1 1997/09/05 19:41:52 jim Exp $";


/*
 * revision history:
 *
 * $Log: dpmta_error.c,v $
 * Revision 1.1  1997/09/05 19:41:52  jim
 * Original distribution.
 *
 * Revision 2.3  1997/05/05  16:17:06  wrankin
 * fixed bad initialization of fmax
 *
 * Revision 2.2  1997/03/18  21:24:21  wrankin
 * fixed error computation code
 *
 * Revision 2.1  1996/10/18  17:04:33  wrankin
 * improved PVM message passing structures.
 * made current and cleaned up T3D library code
 * added additional test codes
 * implement new heirarchical make structure
 *
 *
*/

/* include files */
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "dpmta.h"


main( int argc, char *argv[] )

{
   int nparts;               /* number of particles */

   int i, j, k;              /* loop counters */

   double pdiff, pavg, pmax; /* potential difference, average and max */
   double fdiff, favg, fmax; /* force diference, average and max */
   double prms, frms;        /* rms error */
   double pscale, fscale;

   double dx, dy, dz;
   double scale;

   FILE *fpin, *fpout;

   PmtaPartInfoPtr flist1, flist2;



   /* check and load command line arguments */
   if (argc != 4) {
      fprintf(stderr,
         "%s <#parts> <direct-flist> <dpmta-flist>\n", argv[0]);
      exit(-1);
   }

   nparts = atoi(argv[1]);

   /* allocate particle arrays */
   flist1 = (PmtaPartInfoPtr)malloc(nparts*sizeof(PmtaPartInfo));
   flist2 = (PmtaPartInfoPtr)malloc(nparts*sizeof(PmtaPartInfo));

   /* read in the particle data */
   fpin = fopen(argv[2],"r");
   for ( i=0; i<nparts; i++ ) {
      fscanf(fpin, "%lf %lf %lf %lf ",
         &(flist1[i].f.x), &(flist1[i].f.y), &(flist1[i].f.z), &(flist1[i].v));
/*       fprintf(stderr, "%lf %lf %lf %lf \n", */
/*          flist1[i].f.x, flist1[i].f.y, flist1[i].f.z, flist1[i].v); */
/*       fflush(stderr); */
   }
   fclose(fpin);

   fpin = fopen(argv[3],"r");
   for ( i=0; i<nparts; i++ ) {
      fscanf(fpin, "%lf %lf %lf %lf ",
         &(flist2[i].f.x), &(flist2[i].f.y), &(flist2[i].f.z), &(flist2[i].v));
/*       fprintf(stderr, "%lf %lf %lf %lf \n", */
/*          flist2[i].f.x, flist2[i].f.y, flist2[i].f.z, flist2[i].v); */
/*       fflush(stderr); */
   }
   fclose(fpin);


   /*
   *  print results
   */

   pscale = pavg = prms = pmax = 0.0;
   fscale = favg = frms = fmax = 0.0;

   for ( i=0; i<nparts; i++ ) {
      pdiff = fabs(flist1[i].v - flist2[i].v);
      pscale += fabs(flist2[i].v);
      pavg += pdiff;
      prms += pdiff * pdiff;
      if ( pmax < pdiff ) pmax = pdiff;

      dx = flist1[i].f.x - flist2[i].f.x;
      dy = flist1[i].f.y - flist2[i].f.y;
      dz = flist1[i].f.z - flist2[i].f.z;
      fdiff = sqrt(dx*dx + dy*dy + dz*dz);
      fscale += sqrt(flist2[i].f.x * flist2[i].f.x +
		    flist2[i].f.y * flist2[i].f.y +
		    flist2[i].f.z * flist2[i].f.z );
      favg += fdiff;
      frms += fdiff * fdiff;
      if ( fmax < fdiff ) fmax = fdiff;
   }

   scale = 1.0 / (double)nparts;
   pscale *= scale;
   fscale *= scale;
   pavg *= scale;
   favg *= scale;

   prms *= scale;
   frms *= scale;
   prms = sqrt(prms);   
   frms = sqrt(frms);   


   printf("Avg potential = %lg\n",pscale);
   printf("Avg potential error = %lg (normalized = %lg)\n",pavg, pavg/pscale);
   printf("RMS potential error = %lg (normalized = %lg)\n",prms, prms/pscale);
   printf("Max potential error = %lg (normalized = %lg)\n",pmax, pmax/pscale);

   printf("Avg force |F| = %lg\n",fscale);
   printf("Avg force error |F| = %lg (normalized = %lg)\n",favg, favg/pscale);
   printf("RMS force error |F| = %lg (normalized = %lg)\n",frms, frms/pscale);
   printf("Max force error |F| = %lg (normalized = %lg)\n",fmax, fmax/pscale);

   exit(0);
}

