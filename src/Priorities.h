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

#ifndef PRIORITIES_H
#define PRIORITIES_H

class Priorities
{
public:

static int numBits;

enum PriorityLevels // lower numbers are higher priority
{
	urgent = 0,		// reductions and broadcasts
	high = 1,		// get data messages out fast
	nonlocal = 64,		// computation with off-node data returns
	base = 128,		// default computation priority
	low = 255		// low priority - for big Computes
	/*
	urgent = 0,		// reductions and broadcasts
	high = 0,		// get data messages out fast
	nonlocal = 0,		// computation with off-node data returns
	base = 0,		// default computation priority
	low = 0			// low priority
	*/
};

};

#endif


/***************************************************************************
 * RCS INFORMATION:
 *
 *	$RCSfile $
 *	$Author $	$Locker:  $		$State: Exp $
 *	$Revision: 1.2 $	$Date: 1997/04/06 22:45:11 $
 *
 ***************************************************************************
 * REVISION HISTORY:
 *
 * $Log: Priorities.h,v $
 * Revision 1.2  1997/04/06 22:45:11  ari
 * Add priorities to messages.  Mods to help proxies without computes.
 * Added quick enhancement to end of list insertion of ResizeArray(s)
 *
 * Revision 1.1  1997/04/03 23:22:23  jim
 * Added basic priority() method to Compute.  Only distinguishes between
 * local and nonlocal computations for now.
 *
 *
 ***************************************************************************/
