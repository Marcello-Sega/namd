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

enum PriorityLevels // lower numbers are higher priority
{
	urgent = 0,		// reductions and broadcasts
	nonlocal = 128,		// computation with off-node data returns
	base = 255,		// default computation priority
	low = 512		// low priority
};

};

#endif


/***************************************************************************
 * RCS INFORMATION:
 *
 *	$RCSfile $
 *	$Author $	$Locker:  $		$State: Exp $
 *	$Revision: 1.1 $	$Date: 1997/04/03 23:22:23 $
 *
 ***************************************************************************
 * REVISION HISTORY:
 *
 * $Log: Priorities.h,v $
 * Revision 1.1  1997/04/03 23:22:23  jim
 * Added basic priority() method to Compute.  Only distinguishes between
 * local and nonlocal computations for now.
 *
 *
 ***************************************************************************/
