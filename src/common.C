/**
***  Copyright (c) 1995, 1996, 1997, 1998, 1999, 2000 by
***  The Board of Trustees of the University of Illinois.
***  All rights reserved.
**/

/*
   global functions as declared in common.h
*/

#include "charm++.h"
#include "converse.h"

#include <sys/stat.h>
#include <ctype.h>

#include "common.h"
#include "Communicate.h"
#include "SimParameters.h"
#include "Node.h"


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
   if ( die_hard ) CmiAbort("NAMD ABORTING DUE TO HARD NAMD_quit().\n");
   else CmiAbort("NAMD ABORTING DUE TO SOFT NAMD_quit().\n");
}


// signal all nodes, it's time to quit
void NAMD_die(const char *err_msg)

{
   CkPrintf("FATAL ERROR: %s\n",err_msg);
   CmiAbort(err_msg);
}


// same as write, only does error checking internally
void NAMD_write(int fd, const void *buf, size_t count)

{
   if ( write(fd,buf,count) < 0 ) {
     // NAMD_die("NAMD_write - write to file descriptor failed.");
     NAMD_die(strerror(errno));
   }
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

extern "C" int NAMD_compare_ints(const void *i, const void *j)

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

int *NAMD_bsearch(int *key, int *base, int n, int size, CCMPFN *cmpfn)
{
	if (cmpfn)
		return (int *) bsearch((void *) key,(void *) base,n,size,cmpfn);
	else
	{
		if (*key < base[0] || *key > base[n - 1])
			return NULL;
		if (*key == *base)
			return base;
		int done = 0;
		int top = n;
		int bottom = 0;
		int step = n / 2;
		int i = 0;
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

