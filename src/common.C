/***************************************************************************/
/*    (C) Copyright 1995,1996,1997 The Board of Trustees of the            */
/*                          University of Illinois                         */
/*                           All Rights Reserved                           */
/***************************************************************************/
/***************************************************************************
 * DESCRIPTION:
 *
 * global functions as declared in common.h
 *
 ***************************************************************************/

static char ident[] = "@(#)$Header: /home/cvs/namd/cvsroot/namd2/src/common.C,v 1.1010 1998/02/27 00:14:30 milind Exp $";

#include "chare.h"
#include "ckdefs.h"
#include "c++interface.h"
#include <sys/stat.h>
#include <ctype.h>

#include "common.h"
#include "Communicate.h"
#include "SimParameters.h"
#include "Node.h"

#ifdef GLOBALS
void * operator new (size_t, void *p) { return p; }
#else
void * ::operator new (size_t, void *p) { return p; }
#endif // GLOBALS

#if 0
#ifdef GLOBALS
void * operator new (size_t t) { return CmiAlloc (t); }
void   operator delete (void *p) { if ( p ) CmiFree (p); }
#else
void * ::operator new (size_t t) { return CmiAlloc (t); }
void   ::operator delete (void *p) { if ( p ) CmiFree (p); }
#endif // GLOBALS
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
   cout << "Something called NAMD_quit: die_hard = " << die_hard << endl;
   abort();
}


// signal all nodes, it's time to quit
void NAMD_die(char *err_msg)

{
    CPrintf("%s\n",err_msg);
    abort(); CharmExit();
}


//  NAMD_random returns a random number from 0.0 to 1.0.  It is really
//  just used for portability so porting to machines with different
//  random number generators is easy.  This is less necessary as
//  drand48 appears more and more places . . .

BigReal NAMD_random()
{
	static int first=1;	//  Flag indicating first call

	//  If this is the first call, seed it
	if (first)
	{
		int i;		//  Loop counter

		//  Seed the random number generator initially.  Add the node
		//  number to get more randomness
		srand48(Node::Object()->simParameters->randomSeed+CMyPe());

		//  If this is a client node, now pick off the node-1
		//  first random numbers and then use the next random number
		//  to again seed the generator.  This will really insure that
		//  there is no correlation between random numbers on different
		//  nodes.
		if (CMyPe())
		{
			for (i=0; i<(CMyPe()-1); i++)
			{
				drand48();
			}

			srand48((long) (drand48()*1073741824));
		}

		first=0;
	}

	return(drand48());
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

/***************************************************************************
 Fopen(char *Filename, char *mode):  similar to fopen(filename,mode) except
 it checks for compressed file names too.
 For example:  Fopen("config");
   This will first look for the filename "config" (and return "r" file handle
   if it is found).
   Then it will look for "config.Z" (and run "zcat config.Z", returning
   a file handle to the uncompressed data if found).
   Then it will look for "config.gz" (and run "gzip -d -c config.gz", returning
   a file handle to the uncompressed data if found).
 ***************************************************************************/
FILE *Fopen	(const char *filename, const char *mode)
{
  struct stat buf;
  // check if basic filename exists (and not a directory)
  if (!stat(filename,&buf))
	{
	if (!S_ISDIR(buf.st_mode))
		return(fopen(filename,mode));
	}
  // check for a compressed file
  char *realfilename;
  char *command;
  FILE *fout;
  command = (char *)malloc(strlen(filename)+25);
  // check for .Z (unix compress)
  sprintf(command,"zcat %s.Z",filename);
  realfilename = command+5;
  iout << "Command = " << command << "\n" << endi;
  iout << "Filename.Z = " << realfilename << "\n" << endi;
  if (!stat(realfilename,&buf))
	{
	if (!S_ISDIR(buf.st_mode))
		{
		fout = popen(command,mode);
		// on HP-UX, the first character(s) out of pipe may be
		// garbage!  (Argh!)
		int C;
		do
		  {
		  C = fgetc(fout);
		  // iout << "C is " << C << "\n" << endi;
		  if (isalnum(C) || isspace(C))
			{
			ungetc(C,fout);
			C = -1;	// outta loop
			}
		  } while(C != -1);
		free(command);
		return(fout);
		}
	}
  // check for .gz (gzip)
  sprintf(command,"gzip -d -c %s.gz",filename);
  realfilename = command+11;
  iout << "Command = " << command << "\n" << endi;
  iout << "Filename.gz = " << realfilename << "\n" << endi;
  if (!stat(realfilename,&buf))
	{
	if (!S_ISDIR(buf.st_mode))
		{
		fout = popen(command,mode);
		// on HP-UX, the first character(s) out of pipe may be
		// garbage!  (Argh!)
		int C;
		do
		  {
		  C = fgetc(fout);
		  // iout << "C is " << C << "\n" << endi;
		  if (isalnum(C) || isspace(C))
			{
			ungetc(C,fout);
			C = -1;	// outta loop
			}
		  } while(C != -1);
		free(command);
		return(fout);
		}
	}
  free(command);
  return(NULL);
} /* Fopen() */

/***************************************************************************
 Fclose(FILE *fout):  similar to fclose(fout) except it first checks if the
 file handle fout is a named pipe.
 ***************************************************************************/
int	Fclose	(FILE *fout)
{
  int rc;
  rc = pclose(fout);
  if (rc == -1)	// stream not associated with a popen()
    {
    rc = fclose(fout);
    }
  return rc;
} /* Fclose() */

/***************************************************************************
 * RCS INFORMATION:
 *
 *	$RCSfile: common.C,v $
 *	$Author: milind $	$Locker:  $		$State: Exp $
 *	$Revision: 1.1010 $	$Date: 1998/02/27 00:14:30 $
 *
 ***************************************************************************
 * REVISION HISTORY:
 *
 * $Log: common.C,v $
 * Revision 1.1010  1998/02/27 00:14:30  milind
 * Reduced memory requirements further by using CmiAlloc and CmiFree only
 * for messages and not for all new and deletes.
 *
 * Revision 1.1009  1997/12/26 23:11:06  milind
 * Made namd2 to compile, link and run under linux. Merged Templates and src
 * directoriies, and removed separate definition and declaration files for
 * templates.
 *
 * Revision 1.1008  1997/04/10 18:44:34  nealk
 * 1. changed endl to endi on Controller.C
 * 2. identified popen() bug under HP-UX 9.  popen() occasionally (1/3 of the
 * time) includes garbage characters at the head of the stream, and may
 * close the stream prematurely.  I corrected the garbage bug.  Still need
 * to correct for the closing bug.  Ugh.
 *
 * Revision 1.1007  1997/04/07 14:54:39  nealk
 * Changed fclose() to Fclose() (found in common.[Ch]) to use with popen().
 * Also corrected compilation warnings in Set.[Ch].
 *
 * Revision 1.1006  1997/04/04 17:31:44  brunner
 * New charm fixes for CommunicateConverse, and LdbCoordinator data file
 * output, required proxies, and idle time.
 *
 * Revision 1.1005  1997/04/03 19:59:14  nealk
 * 1) New Fopen() which handles .Z and .gz files.
 * 2) localWaters and localNonWaters lists on each patch.
 *
 * Revision 1.1004  1997/03/19 11:54:59  ari
 * Add Broadcast mechanism.
 * Fixed RCS Log entries on files that did not have Log entries.
 * Added some register variables to Molecule and ComputeNonbondedExcl.C
 *
 * Revision 1.1003  1997/02/13 22:27:05  jim
 * Added inital velocity code from NAMD 1.
 * Reading velocity pdb file appears to work.
 * Reading binary velociy file should work but is untested.
 * Random velocites appears to work but differs from NAMD 1.
 *
 * Revision 1.1002  1997/02/13 16:17:22  ari
 * Intermediate debuging commit - working to fix deep bug in migration?
 *
 * Revision 1.1001  1997/02/11 18:52:00  ari
 * Modified with #ifdef DPMTA to safely eliminate DPMTA codes
 * fixed non-buffering of migration msgs
 * Migration works on multiple processors
 *
 * Revision 1.1000  1997/02/06 15:59:33  ari
 * Resetting CVS to merge branches back into the main trunk.
 * We will stick to main trunk development as suggested by CVS manual.
 * We will set up tags to track fixed points of development/release
 * as suggested by CVS manual - all praise the CVS manual.
 *
 * Revision 1.778  1997/01/28 00:31:36  ari
 * internal release uplevel to 1.778
 *
 * Revision 1.777  1997/01/17 19:37:10  ari
 * Internal CVS leveling release.  Start development code work
 * at 1.777.1.1.
 *
 * Revision 1.5  1996/12/19 17:49:34  jim
 * added null-pointer check to ::delete
 *
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
