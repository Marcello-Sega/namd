/* 
*  dpmta_test2b.c - slave PVM process for DPMTA implentation
*     useing the distributed calling scheme.
*
*  w. t. rankin
*
*  Copyright (c) 1994 Duke University
*  All rights reserved
*
*
*/

static char rcsid[]="$Id: dpmta_test2b.c,v 1.3 1997/09/29 23:58:45 jim Exp $";

/*
 * revision history:
 *
 * $Log: dpmta_test2b.c,v $
 * Revision 1.3  1997/09/29 23:58:45  jim
 * Incorporated changes from version 2.6.1 of DPMTA.
 *   - fixes for bad handling of empty/invalid multipoles when
 *     using large processor sets.
 *   - moved functions that provide data mapping to processors.  master
 *     and slave routines now call the same function in dpmta_distmisc.c
 * Also, switched pvmc.h back to pvm3.h.
 *
 * Revision 1.2  1997/09/12 22:56:36  jim
 * Modifications to work with converse pvm.
 *
 * Revision 1.1  1997/09/05 19:42:12  jim
 * Original distribution.
 *
 * Revision 2.3  1997/05/07  21:28:02  chumphre
 * cubectr changed to cellctr and other minor changes related to the interface
 *
 * Revision 2.2  1996/09/24  18:44:15  wrankin
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
 * Revision 2.1  1995/06/13  04:26:24  wrankin
 * No updates.
 * Sources placed under CVS.
 *
 * Revision 1.2  1994/12/06  19:18:44  wrankin
 * fixed generation of random particle positions
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
#include "pvm3.h"
#include "dpmta.h"

main( int argc, char *argv[] )
{

   int mytid;                /* my task id */
   int mastid;               /* master's task id */

   int i, j, k;              /* loop counters */

   int iter;                 /* number of iterations */
   int num_parts;            /* total number of particles */

   PmtaVector v1;       /* length of the cell cube sides */
   PmtaVector v2;       /* length of the cell cube sides */
   PmtaVector v3;       /* length of the cell cube sides */
   PmtaVector cellctr;       /* coordinates of cube center */

   PmtaParticlePtr particlelist; /* array of particle information */
   PmtaPartInfoPtr forcelist;    /* array of resulting force information */
   PmtaPartInfoPtr flist_lj;     /* array of resulting force information */


   /* enroll master in pvm */
   mytid = pvm_mytid();
   mastid = pvm_parent();

   /* get initial message from master */
   pvm_recv(mastid,1);
   pvm_upkint(&iter,1,1);
   pvm_upkint(&num_parts,1,1);
   pvm_upkdouble(&(v1.x),1,1);
   pvm_upkdouble(&(v2.y),1,1);
   pvm_upkdouble(&(v3.z),1,1);
   pvm_upkdouble(&(cellctr.x),3,1);

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
   forcelist = (PmtaPartInfoPtr)malloc(num_parts * sizeof(PmtaPartInfo));
   flist_lj = (PmtaPartInfoPtr)malloc(num_parts * sizeof(PmtaPartInfo));

   /* create all particles */

   srand48((long)mytid);

   for (i=0; i<num_parts; i++) {
      particlelist[i].p.x = (drand48() - 0.5) * v1.x + cellctr.x;
      particlelist[i].p.y = (drand48() - 0.5) * v2.y + cellctr.y;
      particlelist[i].p.z = (drand48() - 0.5) * v3.z + cellctr.z;
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
   }

   /*
   *  send response back to master
   */

   pvm_initsend(PvmDataDefault);
   pvm_pkdouble(&(particlelist[0].p.x),num_parts*6,1);
   pvm_pkdouble(&(forcelist[0].f.x),num_parts*4,1);
   pvm_pkdouble(&(flist_lj[0].f.x),num_parts*4,1);
   pvm_send(mastid,2);

   pvm_exit();
   exit(0);
}
