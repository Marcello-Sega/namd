/*
*  dpmta_slvcleanup.c - procedures to deallocate and clean up data
*     structures after all force calculations are complete.  
*
*  w. t. rankin
*
*  Copyright (c) 1994 Duke University
*  All rights reserved
*
*/

static char rcsid[] = "$Id: dpmta_slvclean.c,v 1.1 1997/09/05 19:41:58 jim Exp $";

/*
 * Revision history:
 *
 * $Log: dpmta_slvclean.c,v $
 * Revision 1.1  1997/09/05 19:41:58  jim
 * Original distribution.
 *
 * Revision 2.3  1996/08/20  17:12:39  wrankin
 * Enhancement to allow any number of processes to be utilized.
 * Formally the number of precesses was limited to a power of 2.
 *
 * Revision 2.2  1995/06/27  14:20:21  wrankin
 * Implementation of relative interaction list algorithm.
 * Interaction list generation code removed from master library.
 * Cell Table restructured as a 2-D table indexed by level/cell
 * All direct computations are now performed in a single pass.
 * General cleanup and removal of old files.
 *
 * Revision 2.1  1995/06/13  04:26:08  wrankin
 * No updates.
 * Sources placed under CVS.
 *
 * Revision 1.5  1995/04/09  22:21:06  wrankin
 * renaming of all global variables to prevent conflicts for
 * port of dpmta library to the T3D.
 *
 * Revision 1.4  1995/01/13  16:16:53  wrankin
 * update comments
 *
 * Revision 1.3  1994/11/24  13:54:24  wrankin
 * updates to use static particle structure allocation
 * in preparation for having multiple slave entry points
 *
 * Revision 1.2  1994/10/26  02:44:12  wrankin
 * replaced stub with routine to clean out cell table of all
 * partcle information.
 *
 * Revision 1.1  1994/10/14  04:59:47  wrankin
 * Initial revision
 *
 *
*/

/* include files */

#include <stdio.h>
#include "dpmta_cell.h"
#include "dpmta_slave.h"


/****************************************************************
*
*  Slave_Cleanup() - cleans up the contents of the 
*     the cell table at the end of the force calculations.
*
*  this procedure cleans up the cell table data structures at the 
*  end of one iteration, after the force data has been returned to
*  the calling process and is no longer needed by the slave process.
*  this is done in preparation for the next iteration.
*
*  cleanup is performed by cycling through all cells and setting the
*  particle counts to zero and the multipole valid flags (if any)
*  to false.
*
*  since we use realloc() on the particle/force arrays at each iteration,
*  there is no need to implicitly free() them after each iteration.
*
*/

Slave_Cleanup()

   {

   int i, j, k;
   int num_cells;


   num_cells = Dpmta_LevelLocate[Dpmta_NumLevels];

   for (i=0; i<num_cells; i++) {

      if ( Dpmta_CellTbl[0][i] != NULL ) {

         Dpmta_CellTbl[0][i]->mvalid = FALSE;
         Dpmta_CellTbl[0][i]->n = 0;

         if (Dpmta_CellTbl[0][i]->mdata != NULL ) {
            Dpmta_CellTbl[0][i]->mdata->lvalid = FALSE;

	    } /* if mdata != NULL */
         } /* if CellTbl[0][i] != NULL */
      } /* for i */
   } /* Slave_Cleanup */

