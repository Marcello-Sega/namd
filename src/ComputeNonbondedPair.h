/***************************************************************************/
/*                                                                         */
/*              (C) Copyright 1996 The Board of Trustees of the            */
/*                          University of Illinois                         */
/*                           All Rights Reserved                           */
/*									   */
/***************************************************************************/

/***************************************************************************
 * DESCRIPTION: ComputeNonbonded.h
 *
 ***************************************************************************/

#ifndef COMPUTENONBONDED_H
#define COMPUTENONBONDED_H

#include "ComputePatchPair.h"

class ComputeNonbonded : public ComputePatchPair {

public:
  ComputeNonbonded(ComputeID c, PatchID pid[]) : ComputePatchPair(c,pid) { ; }
  virtual ~ComputeNonbonded();

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
 *	$Revision: 1.3 $	$Date: 1996/10/30 01:16:32 $
 *
 ***************************************************************************
 * REVISION HISTORY:
 *
 * $Log: ComputeNonbondedPair.h,v $
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

