/***************************************************************************/
/*                                                                         */
/*              (C) Copyright 1995 The Board of Trustees of the            */
/*                          University of Illinois                         */
/*                           All Rights Reserved                           */
/*								   	   */
/***************************************************************************/

/***************************************************************************
 * RCS INFORMATION:
 *
 *	$RCSfile: common.C,v $
 *	$Author: milind $	$Locker:  $		$State: Exp $
 *	$Revision: 1.4 $	$Date: 1996/12/11 00:04:23 $
 *
 ***************************************************************************
 * DESCRIPTION:
 *
 * global functions as declared in common.h
 *
 ***************************************************************************
 * REVISION HISTORY:
 *
 * $Log: common.C,v $
 * Revision 1.4  1996/12/11 00:04:23  milind
 * *** empty log message ***
 *
 * Revision 1.3  1996/11/11 19:54:09  nealk
 * Modified to use InfoStream instead of Inform.
 *
 * Revision 1.2  1996/08/15 20:32:14  brunner
 * Made NamdDIE use CPrintf
 *
 * Revision 1.1  1996/08/06 20:38:38  ari
 * Initial revision
 *
 * Revision 1.21  1996/04/15 19:00:28  billh
 * Added NAMD_stringdup routine
 *
 * Revision 1.20  1996/02/22 23:03:43  jean
 * added NAMD_bsearch
 *
 * Revision 1.19  1995/10/13 16:38:14  hazen
 * Memory allocation error-checking added
 *
 * Revision 1.18  1995/07/31  11:25:38  brunner
 * Added feedback form location to title message.
 *
 * Revision 1.17  95/04/10  13:37:50  13:37:50  nelson (Mark T. Nelson)
 * Fixed bug in seed send to srand on remote nodes
 * 
 * Revision 1.16  95/03/08  14:33:48  14:33:48  nelson (Mark T. Nelson)
 * Added copyright
 * 
 * Revision 1.15  95/02/15  16:09:45  16:09:45  nelson (Mark T. Nelson)
 * Changed the way drand48 is seeded on client nodes
 * 
 * Revision 1.14  95/02/14  14:56:30  14:56:30  nelson (Mark T. Nelson)
 * Fixed NAMD_die and NAMD_quit so that namd will really terminate nicely in the
 * case of a failure
 * 
 * Revision 1.13  95/02/11  08:47:14  08:47:14  nelson (Mark T. Nelson)
 * Changed NAMD_die to first free the Node object
 * 
 * Revision 1.12  95/01/30  17:15:14  17:15:14  nelson (Mark T. Nelson)
 * Added namdDebug to check_messages()
 * 
 * Revision 1.11  95/01/26  15:34:46  15:34:46  nelson (Mark T. Nelson)
 * Added definitions for MTS memory system
 * 
 * Revision 1.10  95/01/12  14:29:41  14:29:41  nelson (Mark T. Nelson)
 * Changed so that random number seed comes from SimParameters object
 * 
 * Revision 1.9  94/10/11  10:58:50  10:58:50  nelson (Mark T. Nelson)
 * Added NAMD_compare_ints() function
 * 
 * Revision 1.8  94/10/07  13:41:27  13:41:27  nelson (Mark T. Nelson)
 * Changed seed function for drand48
 * 
 * Revision 1.7  94/09/03  18:12:28  18:12:28  billh (Bill Humphrey)
 * Change to how die routine prints ABORT message, and recast of return
 * value from time() to unsigned int.
 * 
 * Revision 1.6  94/09/01  21:01:36  21:01:36  nelson (Mark T. Nelson)
 * Added Attila's name to author's list
 * 
 * Revision 1.5  94/08/30  14:00:28  14:00:28  nelson (Mark T. Nelson)
 * added routine NAMD_random
 * 
 * Revision 1.4  94/08/11  12:16:07  12:16:07  nelson (Mark T. Nelson)
 * Changed NAMD_quit to use #defined tag READYTODIETAG
 * 
 * Revision 1.3  94/07/03  00:13:57  00:13:57  billh (Bill Humphrey)
 * Updated due to new Communicate and Message class structures.
 * 
 * Revision 1.2  94/06/24  03:09:24  03:09:24  billh (Bill Humphrey)
 * Added NAMD_title, NAMD_quit, NAMD_check_messages global functions.
 * Removed NAMD_warn; using Inform objects now to report information.
 * 
 ***************************************************************************/
static char ident[] = "@(#)$Header: /home/cvs/namd/cvsroot/namd2/src/common.C,v 1.4 1996/12/11 00:04:23 milind Exp $";

#include "chare.h"
#include "ckdefs.h"
#include "c++interface.h"

#include "common.h"
#include "Communicate.h"
#include "SimParameters.h"

//  Provide alternate new and delete operators if the MTS fast allocator is being
//  used
#ifdef MTS
void * ::operator new (size_t t) { return malloc (t); }
void   ::operator delete (void *p) { free (p); }
#else
void * ::operator new (size_t t) { return CmiAlloc (t); }
void   ::operator delete (void *p) { CmiFree (p); }
#endif

// print out title
void NAMD_title(void)

{
}


// make a duplicate of a string
char *NAMD_stringdup(const char *s) {
  char *rs;

  if(!s)
    return NULL;

  rs = new char[strlen(s) + 1];
  strcpy(rs,s);

  return rs;
}


// check for messages to be printed out
void NAMD_check_messages(void)

{
}


// normal program termination
void NAMD_quit(Bool die_hard)

{
   cout << "Something called NAMD_quit" << endl;
   abort();
}


// signal all nodes, it's time to quit
void NAMD_die(char *err_msg)

{
    CPrintf("%s\n",err_msg);
    CharmExit();
}


/********************************************************************************/
/*										*/
/*			FUNCTION NAMD_compare_ints				*/
/*										*/
/*   INPUTS:									*/
/*	i - first integer to compare						*/
/*	j - second integer to compare						*/
/*										*/
/*	This function just compares to integers and returns a -1 if i<j, a 0    */
/*   if i=j, and a 1 if i>j.  It is used to pass to qsort and bsearch routines  */
/*   so that standard sorting and seraching can be done.			*/
/*										*/
/********************************************************************************/

int NAMD_compare_ints(const void *i, const void *j)

{
	if ( *((int *) i) < *((int *) j) )
	{
		return(-1);
	}
	else if ( *((int *) i) == *((int *) j) )
	{
		return(0);
	}
	else 
	{
		return(1);
	}
}

int *NAMD_bsearch(int *key, int *base, int n,
	int size, int (*cmpfn) (const void *, const void *))
// int *key, *base, n, size;
// int (*cmpfn)();
{
	int done, top, bottom, step, i;

/*	fprintf(stdout,"searching for %d in array of size %d\n",*key,n);
 *	for (i = 0; i < n; i++)
 *		fprintf(stdout,"  %d\n",base[i]); */
	if (cmpfn)
		return (int *) bsearch((void *) key,(void *) base,n,size,cmpfn);
	else
	{
		if (*key < base[0] || *key > base[n - 1])
			return NULL;
		if (*key == *base)
			return base;
		done = 0;
		top = n;
		bottom = 0;
		step = n / 2;
		while (!done && step)
		{
			i = bottom + step;
			if (base[i] < *key)
				bottom = i;
			else if (base[i] > *key)
				top = i;
			else
				done = 1;
			step = (top - bottom) / 2;
		}
		return done ? base + i : NULL;
	}
}
