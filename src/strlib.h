//-*-c++-*-
/***************************************************************************/
/*                                                                         */
/*              (C) Copyright 1995 The Board of Trustees of the            */
/*                          University of Illinois                         */
/*                           All Rights Reserved                           */
/*								   	   */
/***************************************************************************/
/***************************************************************************
 * DESCRIPTION:
 *	strlib contains a number of useful routines for doing file I/O
 * and some basic string manipulation.  These routines are used for 
 * reading in the parameter and .psf files
 *
 ***************************************************************************/


#ifndef STRLIB_H

#define STRLIB_H

#include <stdio.h>
#include <ctype.h>
#include <stdlib.h>
#include <string.h>
#include <strings.h>
#include "common.h"

void	NAMD_truncate(char *);		//  Remove trailing spaces from
					//  a string
int	NAMD_read_line(FILE *, char *); //  Read in a line from a file
int	NAMD_find_word(char *, char *); //  Check for given word in a
					//  string
int	NAMD_blank_string(char *);	//  Check to see if a string
					//  is blank
void	NAMD_find_first_word(char *, char *);
					//  Find the first word in a string
int	NAMD_read_int(FILE *, char *);  //  Read an integer from a file
void	NAMD_pad(char *, int);		//  Pad a string with leading spaces
void	NAMD_remove_comment(char *);	//  Remove comments at the end of
					//  a line demarked by !

//  Add definitions for missing library routines in AIX
#ifdef AIX
int strcasecmp(const char s[], const char t[]);
int strncasecmp(const char s[], const char t[], int n);
#endif

#endif
/***************************************************************************
 * RCS INFORMATION:
 *
 *	$RCSfile: strlib.h,v $
 *	$Author: ari $	$Locker:  $		$State: Exp $
 *	$Revision: 1.1001 $	$Date: 1997/03/19 11:55:05 $
 *
 ***************************************************************************
 * REVISION HISTORY:
 *
 * $Log: strlib.h,v $
 * Revision 1.1001  1997/03/19 11:55:05  ari
 * Add Broadcast mechanism.
 * Fixed RCS Log entries on files that did not have Log entries.
 * Added some register variables to Molecule and ComputeNonbondedExcl.C
 *
 * Revision 1.1000  1997/02/06 15:59:39  ari
 * Resetting CVS to merge branches back into the main trunk.
 * We will stick to main trunk development as suggested by CVS manual.
 * We will set up tags to track fixed points of development/release
 * as suggested by CVS manual - all praise the CVS manual.
 *
 * Revision 1.778  1997/01/28 00:31:40  ari
 * internal release uplevel to 1.778
 *
 * Revision 1.777.2.1  1997/01/27 22:45:47  ari
 * Basic Atom Migration Code added.
 * Added correct magic first line to .h files for xemacs to go to C++ mode.
 * Compiles and runs without migration turned on.
 *
 * Revision 1.777  1997/01/17 19:37:16  ari
 * Internal CVS leveling release.  Start development code work
 * at 1.777.1.1.
 *
 * Revision 1.1  1996/08/06 20:38:38  ari
 * Initial revision
 *
 * Revision 1.5  1995/04/10 11:30:38  nelson
 * Added strcasecmp and strncasecmp for AIX
 *
 * Revision 1.4  95/03/08  14:46:25  14:46:25  nelson (Mark T. Nelson)
 * Added copyright
 * 
 * Revision 1.3  95/01/26  14:14:27  14:14:27  nelson (Mark T. Nelson)
 * Added function NAMD_remove_comment for charmm22 parameters
 * 
 * Revision 1.2  94/09/30  09:09:47  09:09:47  nelson (Mark T. Nelson)
 * added NAMD_pad function
 * 
 * Revision 1.1  94/06/22  15:05:57  15:05:57  nelson (Mark T. Nelson)
 * Initial revision
 * 
 ***************************************************************************/
