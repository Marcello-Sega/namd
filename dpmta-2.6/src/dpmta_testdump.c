
/* 
*  dpmta_testdump.c - dump formatted data files
*
*  w. t. rankin
*
*  Copyright (c) 1996 Duke University
*  All rights reserved
*
*/

static char rcsid[]="$Id: dpmta_testdump.c,v 1.1 1997/09/05 19:42:15 jim Exp $";

/*
 * revision history:
 *
 * $Log: dpmta_testdump.c,v $
 * Revision 1.1  1997/09/05 19:42:15  jim
 * Original distribution.
 *
 * Revision 2.4  1997/05/07  21:28:11  chumphre
 * cubectr changed to cellctr and other minor changes related to the interface
 *
 * Revision 2.3  1997/05/07  19:31:54  chumphre
 * implement a uniform interface for rectangular and ||-piped cells
 *
 * Revision 2.2  1997/05/07  18:59:49  chumphre
 * major changes to implement ||-piped cells
 *
 * Revision 2.1  1996/10/28  23:02:03  wrankin
 * additions to test routines to provide processing paramenters as
 *   part of output file.
 * dpmta_direct will now perform macroscopic computations, reading
 *   in processing parameters from particle position file.
 *
 *
*/

/* include files */
#include <stdio.h>
#include <stdlib.h>
#include "dpmta.h"


/****************************************************************
*
*  dump resulting position and force vectors
*/

void Dump_Results(
   int nparts,
   PmtaInitData *id,
   PmtaParticlePtr plist,
   PmtaPartInfoPtr flist,
   PmtaPartInfoPtr flist_lj) {

   int i;
   FILE *fout;

   fout = fopen("master-plist","w");

#ifndef PIPED
   fprintf(fout, "CELL_LENGTH: %lg %lg %lg\n",
	   id->v1.x, id->v2.y, id->v3.z);
#else
   fprintf(fout, "PVector1: %lg %lg %lg\n",
           id->v1.x, id->v1.y, id->v1.z);
   fprintf(fout, "PVector2: %lg %lg %lg\n",
           id->v2.x, id->v2.y, id->v2.z);
   fprintf(fout, "PVector3: %lg %lg %lg\n",
           id->v3.x, id->v3.y, id->v3.z);
#endif
   fprintf(fout, "CELL_CENTER: %lg %lg %lg\n",
	   id->cellctr.x, id->cellctr.y, id->cellctr.z);
   fprintf(fout, "PBC/MACRO: %d %d\n", id->pbc, id->kterm);
   fprintf(fout, "THETA: %lg\n", id->theta);
   fprintf(fout, "NPARTS: %d\n", nparts);

   for ( i=0; i<nparts; i++ ) {
#ifdef COMP_LJ
      fprintf(fout,"%20.16lf %20.16lf %20.16lf %20.16lf %20.16lf %20.16lf\n",
           plist[i].p.x, plist[i].p.y, plist[i].p.z, plist[i].q,
           plist[i].a, plist[i].b);
#else
      fprintf(fout,"%20.16lf %20.16lf %20.16lf %20.16lf\n",
           plist[i].p.x, plist[i].p.y, plist[i].p.z, plist[i].q);
#endif
   } /* for i */
   fclose(fout);


   fout = fopen("master-flist","w");
   for ( i=0; i<nparts; i++ ) {
      fprintf(fout,"%20.16lf %20.16lf %20.16lf %20.16lf\n",
           flist[i].f.x, flist[i].f.y, flist[i].f.z, flist[i].v);
   } /* for i */
   fclose(fout);


#ifdef COMP_LJ
   if ( flist_lj != (PmtaPartInfoPtr)NULL ) {
      fout = fopen("master-lj-flist","w");
      for ( i=0; i<nparts; i++ ) {
         fprintf(fout,"%20.16lf %20.16lf %20.16lf %20.16lf\n",
            flist_lj[i].f.x, flist_lj[i].f.y, flist_lj[i].f.z, flist_lj[i].v);
      } /* for i */
      fclose(fout);
   } /* if flist_lj */
#endif

} /* Dump_Results */

