/***************************************************************************/
/*                                                                         */
/*              (C) Copyright 1996 The Board of Trustees of the            */
/*                          University of Illinois                         */
/*                           All Rights Reserved                           */
/*									   */
/***************************************************************************/

/***************************************************************************
 * DESCRIPTION: ComputeNonbondedPair.h
 *
 ***************************************************************************/

#ifndef COMPUTENONBONDEDPAIR_H
#define COMPUTENONBONDEDPAIR_H

#include "ComputePatchPair.h"

class ComputeNonbondedPair : public ComputePatchPair {

public:
  ComputeNonbondedPair(ComputeID c, PatchID pid[]) : ComputePatchPair(c,pid) { ; }
  // virtual ~ComputeNonbondedPair();

protected :
  // virtual void mapReady() { ComputePatchPair::mapReady(); }
  virtual void doForce(Position* p[2], Force* f[2], AtomProperties* a[2]);

};

#endif
/***************************************************************************
 * RCS INFORMATION:
 *
 *	$RCSfile: ComputeNonbondedPair.h,v $
 *	$Author: jim $	$Locker:  $		$State: Exp $
 *	$Revision: 1.5 $	$Date: 1996/11/05 22:40:05 $
 *
 ***************************************************************************
 * REVISION HISTORY:
 *
 * $Log: ComputeNonbondedPair.h,v $
 * Revision 1.5  1996/11/05 22:40:05  jim
 * commented out undefined virtual destructor
 *
 * Revision 1.4  1996/10/31 21:43:29  jim
 * First incarnation as ...Pair
 *
 * Revision 1.3  1996/10/30 01:16:32  jim
 * added AtomProperties structure in Patch plus boxes, passing, etc.
 *
 * Revision 1.2  1996/10/30 00:16:16  jim
 * Removed PositionArray usage.
 *
 * Revision 1.1  1996/10/29 23:55:54  jim
 * Initial revision
 *
 *
 ***************************************************************************/

