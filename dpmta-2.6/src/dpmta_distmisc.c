/* 
*  dpmta_distmisc.c - misc routines to perform mapping between data and
*    distributed processed id.
*
*  these routines are used by both the master and slave processes.  thus
*  they cannot access any global data structures.  most of these
*  functions were taken from dpmta_slvmisc.c
*
*
*  w. t. rankin
*
*  Copyright (c) 1997 Duke University
*  All rights reserved
*
*/

static char rcsid[] = "$Id: dpmta_distmisc.c,v 1.1 1997/09/29 23:58:36 jim Exp $";

/*
 * revision history:
 *
 * $Log: dpmta_distmisc.c,v $
 * Revision 1.1  1997/09/29 23:58:36  jim
 * Incorporated changes from version 2.6.1 of DPMTA.
 *   - fixes for bad handling of empty/invalid multipoles when
 *     using large processor sets.
 *   - moved functions that provide data mapping to processors.  master
 *     and slave routines now call the same function in dpmta_distmisc.c
 * Also, switched pvmc.h back to pvm3.h.
 *
 * Revision 2.1  1997/09/29 20:24:52  wrankin
 * fixed problem with invalid (empty) multipoles during upward pass.
 * cell indexing by processor was inconsistant between master/slave.
 *
 *
 */

/* include files */
#include <stdio.h>


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
*  returns the pid for the slave that owns a given cell
*/

int getslvpid(int np, int level, int cell)
{
   return ( (np * cell) / (0x01 <<(3*level)) );
}


/*****************************************************************
*
*  returns the first cell that a specific pid owns for a given
*  level.
*/

int getscell( int np, int pid, int level )
{
   int sc, ec;

   sc = (pid * (0x1 << (3*level)) + (np-1)) / np ;
   ec = ((pid+1) * (0x1 << (3*level)) + (np-1)) / np - 1 ;

   if ( sc > ec )
      return(-1);
   else
      return(sc);
}


/*****************************************************************
*
*  returns the last cell that a specific pid owns for a given
*  level.
*/

int getecell( int np, int pid, int level )
{
   int sc, ec;

   sc = (pid * (0x1 << (3*level)) + (np-1)) / np ;
   ec = ((pid+1) * (0x1 << (3*level)) + (np-1)) / np - 1 ;

   if ( sc > ec )
      return(-1);
   else
      return(ec);
}
