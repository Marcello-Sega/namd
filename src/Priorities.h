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

// use for nonlocal communications
static int comm_urgent;		// reductions and broadcasts
static int comm_high;		// data and results
static int comm_default;	// whatever
static int comm_low;		// collections

// use for local enqueueing
static int comp_default;	// default in Compute base class
static int comp_nonlocal_base;	// small patch count, some off-node
static int comp_nonlocal_range;	//  range for above
static int comp_local_large;	// large patch count, all local
static int comp_local_base;	// small patch count, all local
static int comp_local_range;	//  range for above
static int comp_synchronizing;	// causes multi-node synchronization

};

#endif


/***************************************************************************
 * RCS INFORMATION:
 *
 *	$RCSfile $
 *	$Author $	$Locker:  $		$State: Exp $
 *	$Revision: 1.3 $	$Date: 1997/08/26 16:26:16 $
 *
 ***************************************************************************
 * REVISION HISTORY:
 *
 * $Log: Priorities.h,v $
 * Revision 1.3  1997/08/26 16:26:16  jim
 * Revamped prioritites for petter performance and easier changes.
 *
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
