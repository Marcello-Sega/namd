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

// all methods are inlined so far...

// some suggested usages:
// Send a warning message to the console:
//    iout << iWARN << blah << endc;
// Send running information to the information monitor (CPrintf)
//    iout << iINFO << blah << endi;
// Send an error message to the information monitor
//    iout << iERROR << blah << endi;
// Set your own information to level 4
//    iout << blah << level(4) << blah << endi;
// Display all information that is at level 6 or higher
//    iout << blah << endi(6);

// Debug messages should use DebugM()!  not iout.


/***************************************************************************
 * RCS INFORMATION:
 *
 *	$RCSfile $
 *	$Author $	$Locker:  $		$State: Exp $
 *	$Revision: 1.1001 $	$Date: 1997/03/19 11:54:18 $
 *
 ***************************************************************************
 * REVISION HISTORY:
 *
 * $Log: InfoStream.C,v $
 * Revision 1.1001  1997/03/19 11:54:18  ari
 * Add Broadcast mechanism.
 * Fixed RCS Log entries on files that did not have Log entries.
 * Added some register variables to Molecule and ComputeNonbondedExcl.C
 *
 ***************************************************************************/
