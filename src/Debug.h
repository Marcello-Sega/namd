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
#include "ckdefs.h"	// for CPrintf

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
	    if (level >= STDERR_LEVEL)	CPrintf("ERROR (%d): ",level); \
	    else if (level > 0) CPrintf("Debug (%d): ",level); \
	    char debugBuf[256]; \
	    ostrstream dout(debugBuf,sizeof(debugBuf)); \
	    dout << format << ends; \
	    CPrintf("%s %d: %s",__FILE__,__LINE__,debugBuf); \
	  } \
	}

 #else
  /* make a void function. */
  /* parameters with side effects will be removed! */
  #define Debug(x) ;
  #define DebugM(x,y)	;

 #endif /* DEBUGM */

#endif /* DEBUG_H */

