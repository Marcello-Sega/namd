/* 
*  dpmta_slvmisc.c - misc routines to perform slave indexing and other
*                   identification tasks
*
*  w. t. rankin
*  w. elliott
*
*  Copyright (c) 1994 Duke University
*  All rights reserved
*
*/

static char rcsid[] = "$Id: dpmta_slvmisc.c,v 1.1 1997/09/05 19:42:02 jim Exp $";

/*
 * revision history:
 *
 * $Log: dpmta_slvmisc.c,v $
 * Revision 1.1  1997/09/05 19:42:02  jim
 * Original distribution.
 *
 * Revision 2.9  1996/08/20  17:12:49  wrankin
 * Enhancement to allow any number of processes to be utilized.
 * Formally the number of precesses was limited to a power of 2.
 *
 * Revision 2.8  1996/02/29  21:13:46  wrankin
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
 * Revision 2.7  1995/10/01  21:46:14  wrankin
 * implemented LJ multipole calculations - compile time option only
 *
 * Revision 2.6  1995/09/18  19:02:45  wrankin
 * updates to move the scaling of the simulation box to the slave
 * process.
 *
 * Revision 2.5  1995/06/27  14:20:27  wrankin
 * Implementation of relative interaction list algorithm.
 * Interaction list generation code removed from master library.
 * Cell Table restructured as a 2-D table indexed by level/cell
 * All direct computations are now performed in a single pass.
 * General cleanup and removal of old files.
 *
 * Revision 2.4  1995/04/09  22:21:06  wrankin
 * renaming of all global variables to prevent conflicts for
 * port of dpmta library to the T3D.
 *
 * Revision 2.3  1995/02/24  21:22:38  wrankin
 * fixed declaration of getnpids so it correctly returns integer.
 *
 * Revision 2.2  1994/10/19  00:00:46  wrankin
 * added copyright notice
 *
 * Revision 2.1  1994/06/02  12:52:27  wrankin
 * No change.  Update version number for Release 2
 *
 * Revision 1.4  1994/03/18  17:47:21  wrankin
 * renamed getpid() call so no conflict with unix system call
 *
 * Revision 1.3  1994/03/16  14:23:46  wrankin
 * added <stdio.h> header file
 *
 * Revision 1.2  1994/03/13  04:06:24  wrankin
 * update to use global variables
 *
 * Revision 1.1  1994/03/10  16:19:35  wrankin
 * Initial revision
 *
 */

/* include files */
#include <stdio.h>
#include "dpmta_cell.h"
#include "dpmta_slave.h"


/****************************************************************
*
*  returns the parent cell id for any given cell
*/

int getparent( int cellid )
{
   return (cellid >> 3);
} /* get parent */


/****************************************************************
*
*  returns the lowest indexed child cell id for any given cell
*/

int getfirstchild( int cellid )
{
   return (cellid << 3);
} /* get first child */



/****************************************************************
*
*  returns the parent process id for a slave at a given level
*  the returned value may be either the current Pid, indicating
*  that the current slave owns the parent cell at the next level
*  up the tree, or the function will return the Pid of the process
*  that owns the parent of the cells owned by the current pid.
*
*  clear as mud, right?
*/

int parentpid( int level )
{
   int incr;

/*    if (incr = Dpmta_Nproc/Dpmta_Power8[level]) */
/*       return ((Dpmta_Pid/incr) * incr); */
/*    else */
/*       return (Dpmta_Pid); */

   if (level == 0)
      return 0;

   /* number of procs/cell at the parent level */
   incr = Dpmta_Nproc/Dpmta_Power8[level-1];
   if (incr <= 1)
      return (Dpmta_Pid);
   else
      return ((Dpmta_Pid/incr) * incr);

}


/****************************************************************
*
*  returns the pid for the slave that owns a given cell
*/

int getslvpid(int level, int cell)
{
   return ( (Dpmta_Nproc * cell) / (0x01 <<(3*level)) );
}


/****************************************************************
*
*  returns the number of pids(slaves) for a given level
*/

int getnpids(int level)
{
   int npids;

   npids = Dpmta_Power8[level];
   if (npids <= Dpmta_Nproc)
      return (npids);
   else
      return (Dpmta_Nproc);
}


