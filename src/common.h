/**
***  Copyright (c) 1995, 1996, 1997, 1998, 1999, 2000 by
***  The Board of Trustees of the University of Illinois.
***  All rights reserved.
**/

/*
   common definitions for namd.
*/

#ifndef COMMON_DEF_H
#define COMMON_DEF_H

#include <unistd.h>
#include <stdio.h>
#include <stdlib.h>
#include <iostream.h>
#include <errno.h>
#include <limits.h>
#include "InfoStream.h"

#if ( INT_MAX == 2147483647 )
typedef	int	int32;
#else
typedef	short	int32;
#endif

#if defined(PLACEMENT_NEW)
void * ::operator new (size_t, void *p) { return p; }
#elif defined(PLACEMENT_NEW_GLOBAL)
void * operator new (size_t, void *p) { return p; }
#endif

#ifndef DEFPRIO
#define DEFPRIO (1 * 64)
#endif

#define COLOUMB 332.0636
#define BOLTZMAN 0.001987191
#define TIMEFACTOR 48.88821
#define PRESSUREFACTOR 6.95E4
#define PDBVELFACTOR 20.45482706
#define PDBVELINVFACTOR (1.0/PDBVELFACTOR)
#define PNPERKCALMOL 69.479

#ifndef PI
#define PI	3.141592653589793
#endif

#ifndef TWOPI
#define TWOPI	2.0 * PI
#endif

#ifndef ONE
#define ONE	1.000000000000000
#endif

#ifndef ZERO
#define ZERO	0.000000000000000
#endif

#ifndef SMALLRAD
#define SMALLRAD      0.0005
#endif

#ifndef SMALLRAD2
#define SMALLRAD2     SMALLRAD*SMALLRAD
#endif

/* Define the size for Real and BigReal.  Real is usually mapped to float */
/* and BigReal to double.  To get BigReal mapped to float, use the 	  */
/* -DSHORTREALS compile time option					  */
typedef float	Real;

#ifdef SHORTREALS
typedef float	BigReal;
#else
typedef double  BigReal;
#endif

#ifndef FALSE
#define FALSE 0
#define TRUE 1
#endif

#ifndef NO
#define NO 0
#define YES 1
#endif

#ifndef STRINGNULL
#define STRINGNULL '\0'
#endif

#define MAX_NEIGHBORS 27

#define NAMD_DEBUG	FALSE

typedef int Bool;

// forward declarations of global variables
class Inform;
extern Inform namdInfo;
extern Inform namdWarn;
extern Inform namdErr;
extern Inform namdDebug;

class Communicate;

// global functions
void NAMD_title(void);
void NAMD_check_messages(void);
void NAMD_quit(Bool die_hard=FALSE);
void NAMD_die(const char *);
void NAMD_bug(const char *);
void NAMD_write(int fd, const void *buf, size_t count); // NAMD_die on error
char *NAMD_stringdup(const char *);
extern "C" {
  int NAMD_compare_ints(const void *, const void *);
  typedef int CCMPFN(const void *, const void*);
}
int *NAMD_bsearch(int *, int *, int, int, CCMPFN *cmpfn);
FILE *Fopen(const char *filename, const char *mode);
int  Fclose(FILE *fout);

// message tags
// NOTE!!!  Do NOT use any tags smaller than 100.  Add tags sequentially
// and update the MAXTAGVALUE parameter after adding tag
#define SIMPARAMSTAG	100	//  Tag for SimParameters class
#define STATICPARAMSTAG 101	//  Tag for Parameters class
#define MOLECULETAG	102	//  Tag for Molecule class
#define DISTRIBTAG	103	//  Tag for patch distributions
				//  going from host to client nodes
#define INITIALPOSTAG	104	//  Tag for Initial Position message
				//  from host to client nodes
#define INITIALVELTAG	105	//  Tag for Initial Velocities message

#define NBCOORTAG	106	//  Tag for non-bonded coordinates being
				//  sent from one patch to another
#define NBFORCETAG	107	//  Tag for non-bonded forces being
				//  sent from on patch to another
#define BONDCOORTAG	108	//  Tag for bonded coordinates being sent
				//  from one patch to another
#define LOCALCALCTAG	109	//  Tag for the self message indicating
				//  that local calculations should be
				//  performed
#define INTEGRATETAG	110	//  Tag for the self message that triggers
				//  integration
#define MULTINBCOORTAG	111	//  Tag for nb coordinates to multiple patch
				//  message.

#define OUTPUTENERGY    112     // Tag for energy output messages,
                                //  collected from nodes
#define OUTPUTCOOR      113     // Tag for coordinate output messages,
                                //  collected from nodes
#define OUTPUTVEL       114     // Tag for velocity output messages,
                                // collected from nodes 
#define ATOMREASSIGNTAG 115	// Tag for atom reassignment messages
#define PATCHMOVETAG    116     // Tag for a patch movement message

#define PATCHLOADTAG	117	// Tag for patch load data
#define LSTATINFOTAG    118     // Tag to distribute load stats 
				// among processors
#define LDBBARRIERTAG   119     // Wait at barrier for this message

#define RESCALEVELTAG   120     // Temperature rescaling information

#define ERRTAG		121	//  Tag for namdErr
#define WARNTAG		122	//  Tag for namdWarn
#define INFOTAG		123	//  Tag for namdInfo
#define DEBUGTAG	124	//  Tag for namdDebug
#define READYTODIETAG	125	//  Tag used by NAMD_quit
#define DIENOWTAG	126	//  Tag used by NAMD_die

#define COMTAG		127	//  Tag used to send Center of Mass
#define FULLTAG		128	//  Tag for full electrostatics coordinates
#define FULLFORCETAG	129	//  Tag for full electrostatics forces
#define OUTPUTLONGFORCE	130	//  Tag for force DCD generation
#define OUTPUTSHORTFORCE 131	//  Tag for short force DCD generation
#define LASTTEMPTAG     132	//  Tag for last temperature for temperature
				//  coupling
#define OUTPUTALLFORCE  133	//  Tag for total force DCD generation
#define PREINTEGRATETAG	134	//  Tag to send x,v,f to GlobalIntegrate
#define GLOBALINTTAG	135	//  Tag for GlobalIntegrate messages

#define DPMTATAG	136	//  Tag for DPMTA

#define SMDDATATAG      137     //  Tag for sending SMDData to clients

#define MAXTAGVALUE	137	//  Maximum tag value.  NOTE!!  Always
				//  reset this when new tags are added,
				//  otherwise MessageManager object will
				//  almost certainly cause a core dump!!

/*****************************************************************************/
/*									     */
/*			MACRO DEBUG_MSG, DEBUG_MSG_N			     */
/*									     */
/*     This prints a message to the debug Inform object; it is used to allow */
/*  selective inclusion or exclusion of debugging output from the executable.*/
/*     DEBUG_MSG takes one argument, and just prints the message to the debug*/
/*  Inform object, setting the current message level to 1.		     */
/*     DEBUG_MSG_N takes two arguments: the level of the message, and the    */
/*  message to print.  Otherwise it acts like DEBUG_MSG.		     */
/*									     */
/*****************************************************************************/

#if NAMD_DEBUG

#define DEBUG_MSG(m) { if(namdDebug.on()) { namdDebug << m; } }
#define DEBUG_MSG_N(l,m) { if(namdDebug.on()) { namdDebug.msg_level(l) << m; } }

#else

#define DEBUG_MSG(m)
#define DEBUG_MSG_N(l,m)

#endif

/*************************************************************************************/
/*										     */
/*				MACRO CHECK_DOMAIN				     */
/*										     */
/*	This function checks for domain errors from mathmatical functions.  If one   */
/*  is found, then an error message containing the line number and file that the     */
/*  error occured in is printed.						     */
/*										     */
/*************************************************************************************/

#define CHECK_DOMAIN() \
	{ \
		if (errno == EDOM) \
		{ \
		  iout << iERRORF << "Domain error\n" << endi; \
		  NAMD_die("Domain error"); \
		} \
	}

#define CHECK_DOMAIN_ACOS(x) \
	{ \
	if (errno == EDOM) \
	  { \
	    iout << iERRORF << "Domain error with acos(" << x << ")\n" << endi; \
	    NAMD_die("Domain error"); \
	  } \
	}

#define CHECK_DOMAIN_ATAN(x,y) \
	{ \
	if (errno == EDOM) \
	  { \
	    iout << iERRORF << "Domain error with atan2(" << x << "," << y << ")\n" << endi; \
	    NAMD_die("Domain error"); \
	  } \
	}

#endif

