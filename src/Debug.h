//-*-c++-*-
/***************************************************************************/
/*                                                                         */
/*              (C) Copyright 1996 The Board of Trustees of the            */
/*                          University of Illinois                         */
/*                           All Rights Reserved                           */
/*                                                                         */
/***************************************************************************/

/***************************************************************************
 * DESCRIPTION:
 *
 ***************************************************************************/
#include "InfoStream.h"

#ifndef DEBUG_H
#define DEBUG_H

#ifndef MIN_DEBUG_LEVEL
  #define MIN_DEBUG_LEVEL 0
#endif
#ifndef MAX_DEBUG_LEVEL
  #define MAX_DEBUG_LEVEL 10
#endif
#ifndef STDERR_LEVEL
  /* anything >= this error level goes to stderr */
  #define STDERR_LEVEL 5
#endif

#include <stdio.h>
#include <stdarg.h>
#include <strstream.h>	// for ostrstream
#include "charm++.h"	// for CkPrintf

/*****************************************************************
 *  DebugM(): function to display a debug message.
 *  Messages have different levels.  The low numbers are low severity
 *  while the high numbers are really important.  Very high numbers
 *  are sent to stderr rather than stdout.
 *  The default severity scale is from 0 to 10.
 *     0 = plain message
 *     4 = important message
 *     5 = warning (stderr)
 *     10 = CRASH BANG BOOM error (stderr)
 *  The remaining args are like printf: a format string and some args.
 *  This function can be turned off by compiling without the DEBUGM flag
 *  No parameters to this function should have a side effect!
 *  No functions should be passed as parameters!  (including inline)
 *****************************************************************/
 #ifdef DEBUGM

  #define Debug(x) (x)
  #define DebugM(level,format) \
	{ \
	  if ((level >= MIN_DEBUG_LEVEL) && (level <= MAX_DEBUG_LEVEL)) \
	  { \
	    infostream Dout; \
	    if (level >= STDERR_LEVEL)	Dout << "[ERROR " << level << "] "; \
	    else if (level > 0) Dout << "[Debug " << level << "] "; \
	    Dout << iPE << ' ' << iFILE; \
	    Dout << format << endi; \
	  } \
	}

 #else
  /* make a void function. */
  /* parameters with side effects will be removed! */
  #define Debug(x) ;
  #define DebugM(x,y)	;

 #endif /* DEBUGM */

#endif /* DEBUG_H */


/***************************************************************************
 * RCS INFORMATION:
 *
 *	$RCSfile $
 *	$Author $	$Locker:  $		$State: Exp $
 *	$Revision: 1.1003 $	$Date: 1999/05/11 23:56:31 $
 *
 ***************************************************************************
 * REVISION HISTORY:
 *
 * $Log: Debug.h,v $
 * Revision 1.1003  1999/05/11 23:56:31  brunner
 * Changes for new charm version
 *
 * Revision 1.1002  1998/03/03 23:05:12  brunner
 * Changed include files for new simplified Charm++ include file structure.
 *
 * Revision 1.1001  1997/03/19 11:54:15  ari
 * Add Broadcast mechanism.
 * Fixed RCS Log entries on files that did not have Log entries.
 * Added some register variables to Molecule and ComputeNonbondedExcl.C
 *
 ***************************************************************************/
