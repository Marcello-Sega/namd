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
#ifndef SP2
#include <strings.h>
#endif
#include "common.h"

void	NAMD_truncate(char *);		//  Remove trailing spaces from
					//  a string
int	NAMD_read_line(FILE *, char *); //  Read in a line from a file
int	NAMD_find_word(const char *, const char *); //  Check for given word in a
					//  string
int	NAMD_blank_string(char *);	//  Check to see if a string
					//  is blank
void	NAMD_find_first_word(char *, char *);
					//  Find the first word in a string
int	NAMD_read_int(FILE *, const char *);  //  Read an integer from a file
void	NAMD_pad(char *, size_t);	//  Pad a string with leading spaces
void	NAMD_remove_comment(char *);	//  Remove comments at the end of
					//  a line demarked by !

//  Add definitions for missing library routines in AIX
#ifdef SP2
int strcasecmp(const char s[], const char t[]);
int strncasecmp(const char s[], const char t[], int n);
#endif

#endif

