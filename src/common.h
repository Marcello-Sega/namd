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
 *	$RCSfile: common.h,v $
 *	$Author: nealk $	$Locker:  $		$State: Exp $
 *	$Revision: 1.2 $	$Date: 1996/12/02 17:05:37 $
 *
 ***************************************************************************
 * DESCRIPTION:
 *
 * common definitions for namd.
 ***************************************************************************
 * REVISION HISTORY:
 *
 * $Log: common.h,v $
 * Revision 1.2  1996/12/02 17:05:37  nealk
 * Debugging stuff.
 *
 * Revision 1.1  1996/08/06 20:38:38  ari
 * Initial revision
 *
 * Revision 1.49  1996/04/23 20:41:54  billh
 * Added SMALLRAD and SMALLRAD2 defines, which indicate minimum distances
 * between atoms before an error occurs in a force or energy evaluation.
 *
 * Revision 1.48  1996/04/15 19:00:28  billh
 * Added NAMD_stringdup routine
 *
 * Revision 1.47  1996/02/22 23:03:43  jean
 * added NAMD_bsearch
 *
 * Revision 1.46  1995/10/29 00:59:05  jim
 * Eliminated DONEPREINTEGRATETAG.
 *
 * Revision 1.45  95/10/27  21:37:40  21:37:40  jim (Jim Phillips)
 * Added tags for global integration.
 * 
 * Revision 1.44  95/09/26  15:15:14  15:15:14  nelson (Mark T. Nelson)
 * Added tag for collection of all forces for DCD output
 * 
 * Revision 1.43  95/09/26  13:27:54  13:27:54  nelson (Mark T. Nelson)
 * Added tag for temperature coupling'
 * 
 * Revision 1.42  95/09/21  17:43:58  17:43:58  billh (Bill Humphrey)
 * Added checks to see if things defined such as PI, TWOPI, TRUE, FALSE, etc
 * 
 * Revision 1.41  95/08/30  11:08:18  11:08:18  nelson (Mark T. Nelson)
 * Added MAXTAGVALUE tag
 * 
 * Revision 1.40  95/08/11  14:06:18  14:06:18  nelson (Mark T. Nelson)
 * Added OUTPUTFORCE tag for force DCD file generation
 * 
 * Revision 1.39  95/05/13  10:43:36  10:43:36  nelson (Mark T. Nelson)
 * Added SHORTREALS ifdef
 * 
 * Revision 1.38  95/04/06  14:44:48  14:44:48  nelson (Mark T. Nelson)
 * Removed definition of FIRST_TIMESTEP
 * 
 * Revision 1.37  95/04/06  12:52:06  12:52:06  nelson (Mark T. Nelson)
 * Removed extern class references
 * 
 * Revision 1.36  95/03/30  21:55:32  21:55:32  nelson (Mark T. Nelson)
 * Added messages for full electrostatics
 * 
 * Revision 1.35  95/03/22  11:23:45  11:23:45  nelson (Mark T. Nelson)
 * Added tag for sending center of mass
 * 
 * Revision 1.34  95/03/18  13:11:56  13:11:56  nelson (Mark T. Nelson)
 * Changed TAGS to work with new MessageManager scheme
 * 
 * Revision 1.33  95/03/10  16:28:59  16:28:59  brunner (Robert Brunner)
 * Added MAX_NEIGHBORS=27 for declaring max array sizes and
 * MULTINBCOORTAG tag for multiple coordinate messages.
 * 
 * Revision 1.32  95/03/08  14:35:03  14:35:03  nelson (Mark T. Nelson)
 * Added copyright
 * 
 * Revision 1.31  95/03/06  20:52:37  20:52:37  nelson (Mark T. Nelson)
 * Added include of errno.h for EDOM definition
 * 
 * Revision 1.30  95/02/14  14:56:55  14:56:55  nelson (Mark T. Nelson)
 * Fixed NAMD_die and NAMD_quit so that namd will really terminate nicely in the
 * case of a failure
 * 
 * Revision 1.29  95/02/01  14:15:57  14:15:57  nelson (Mark T. Nelson)
 * Replaced ifdef with if for NAMD_DEBUG defintion
 * 
 * Revision 1.28  95/02/01  11:49:06  11:49:06  billh (Bill Humphrey)
 * Added DEBUG_MSG, DEBUG_MSG_N macros
 * 
 * Revision 1.27  95/01/30  20:09:32  20:09:32  nelson (Mark T. Nelson)
 * Added RESCALEVELTAG
 * 
 * Revision 1.26  95/01/26  15:34:38  15:34:38  nelson (Mark T. Nelson)
 * Added definitions for MTS memory system
 * 
 * Revision 1.25  95/01/19  13:29:04  13:29:04  brunner (Robert Brunner)
 * Added PATCHMOVETAG, which indicates a message for patch migration
 * due to load balancing.
 * 
 * Revision 1.24  94/12/01  14:16:44  14:16:44  brunner (Robert Brunner)
 * Added LSTATINFOTAG and LDBBARRIERTAG for load balancing
 * 
 * Revision 1.23  94/11/28  14:14:01  14:14:01  nelson (Mark T. Nelson)
 * Fixed so that line and file in CHECK_DOMAIN were correct
 * 
 * Revision 1.22  94/11/23  09:32:10  09:32:10  nelson (Mark T. Nelson)
 * Added macro CHECK_DOMAIN
 * 
 * Revision 1.21  94/11/22  13:39:02  13:39:02  nelson (Mark T. Nelson)
 * Added tag for patch load stats
 * 
 * Revision 1.20  94/10/11  10:59:05  10:59:05  nelson (Mark T. Nelson)
 * Added NAMD_compare_ints() function and ATOMREDISTTAG
 * 
 * Revision 1.19  94/10/07  12:19:51  12:19:51  nelson (Mark T. Nelson)
 * Corrected value of boltzman's constant
 * 
 * Revision 1.18  94/10/06  15:44:15  15:44:15  gursoy (Attila Gursoy)
 * removed the SYNCMSG type (no need anymore)
 * 
 * Revision 1.17  94/09/29  18:27:55  18:27:55  gursoy (Attila Gursoy)
 * defined a new message tag for synchronization SYNCMSG
 * 
 * Revision 1.16  94/09/28  18:29:14  18:29:14  gursoy (Attila Gursoy)
 * added output message tags OUTPUTENERGY, OUTPUTCOOR, OUTPUTVEL
 * 
 * Revision 1.15  94/09/22  12:16:10  12:16:10  nelson (Mark T. Nelson)
 * Added COLOUMB and TIMEFACTOR
 * 
 * Revision 1.14  94/09/09  10:20:01  10:20:01  gursoy (Attila Gursoy)
 * added FIRST_TIMESTEP
 * 
 * Revision 1.13  94/09/02  16:36:24  16:36:24  nelson (Mark T. Nelson)
 * moved NAMD_DEBUG from namd.C to common.h
 * 
 * Revision 1.12  94/08/30  13:59:55  13:59:55  nelson (Mark T. Nelson)
 * Added a bunch of tags
 * 
 * Revision 1.11  94/08/09  13:20:15  13:20:15  nelson (Mark T. Nelson)
 * Added DISTRIBTAG
 * 
 * Revision 1.10  94/08/03  21:59:03  21:59:03  nelson (Mark T. Nelson)
 * Added TAG for Molecule class
 * 
 * Revision 1.9  94/08/02  16:29:03  16:29:03  nelson (Mark T. Nelson)
 * Added tags for Parameters class and Inform objcets
 * 
 * Revision 1.8  94/08/01  10:38:20  10:38:20  nelson (Mark T. Nelson)
 * Changed definitions of Bool, TRUE, FALSE, etc. to work better
 * with the communication class
 * 
 * Revision 1.7  94/07/08  16:37:35  16:37:35  nelson (Mark T. Nelson)
 * Added message tag for SimParameters
 * 
 * Revision 1.6  94/07/06  15:21:11  15:21:11  nelson (Mark T. Nelson)
 * Added BigReal definition
 * 
 * Revision 1.5  94/06/29  22:48:37  22:48:37  dalke (Andrew Dalke)
 * Changed Bool from int to an enum
 * 
 * Revision 1.4  94/06/28  17:42:43  17:42:43  billh (Bill Humphrey)
 * Added namdDebug global variable for debugging messages.  Turned
 * on by defining NAMD_DEBUG.
 * 
 * Revision 1.3  94/06/24  03:10:11  03:10:11  billh (Bill Humphrey)
 * Added NAMD_title, NAMD_quit, NAMD_check_messages global functions.
 * Removed NAMD_warn; using Inform objects now to report information.
 * 
 * Revision 1.2  94/06/22  15:07:20  15:07:20  nelson (Mark T. Nelson)
 * Added NAMD_die and NAMD_warn prototypes
 * 
 * Revision 1.1  94/06/20  21:36:16  21:36:16  billh (Bill Humphrey)
 * Initial revision
 * 
 ***************************************************************************/

#ifndef COMMON_DEF_H
#define COMMON_DEF_H

#include <stdio.h>
#include <stdlib.h>
#include <iostream.h>
#include <errno.h>
#include "InfoStream.h"

//  Redefine new and delete if the MTS fast malloc is being used
#ifdef MTS
void * ::operator new (size_t);
void   ::operator delete (void *);
#endif

#define COLOUMB 332.0636
#define BOLTZMAN 0.001987191
#define TIMEFACTOR 48.88821

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
extern Communicate *comm;

// global functions
void NAMD_title(void);
void NAMD_check_messages(void);
void NAMD_quit(Bool die_hard=FALSE);
void NAMD_die(char *);
int NAMD_compare_ints(const void *, const void *);
char *NAMD_stringdup(const char *);
int *NAMD_bsearch(int *, int *, int, int, int (*cmpfn) (const void *, const void
*));

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

#define MAXTAGVALUE	135	//  Maximum tag value.  NOTE!!  Always
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
