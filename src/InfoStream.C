/***************************************************************************/
/*                                                                         */
/*              (C) Copyright 1996 The Board of Trustees of the            */
/*                          University of Illinois                         */
/*                           All Rights Reserved                           */
/*                                                                         */
/***************************************************************************/

/***************************************************************************
 * DESCRIPTION:
 *  Defines a new stream: iout, for "i"nforming consoles.
 ***************************************************************************/

#include "InfoStream.h"

infostream iout;

/* later, endi and endc should become modifiers! */
#undef endi
#undef endc

char * endi(infostream& s)	{ s.endi(); return ""; }
char * endc(infostream& s)	{ s.endc(); return ""; }


// some suggested usages:
// Send a warning message to the console:
//    iout << iWARN << blah << endc;
// Send running information to the information monitor (CPrintf)
//    iout << iINFO << blah << endi;
// Send an error message to the information monitor
//    iout << iERROR << blah << endi;

// Debug messages should use DebugM()!  not iout.

