/***************************************************************************/
/*                                                                         */
/*              (C) Copyright 1996 The Board of Trustees of the            */
/*                          University of Illinois                         */
/*                           All Rights Reserved                           */
/*									   */
/***************************************************************************/

/***************************************************************************
 * DESCRIPTION: ComputeNonbondedSelf.h
 *
 ***************************************************************************/

#ifndef COMPUTENONBONDEDSELF_H
#define COMPUTENONBONDEDSELF_H

#include "ComputePatch.h"
#include "ComputeNonbondedUtil.h"

class ComputeNonbondedSelf : public ComputePatch, ComputeNonbondedUtil {

public:
  ComputeNonbondedSelf(ComputeID c, PatchID pid);
  virtual ~ComputeNonbondedSelf();

protected :
  // virtual void mapReady() { ComputePatch::mapReady(); }
  virtual void doForce(Position* p, Force* f, AtomProperties* a);

  ReductionMgr *reduction;

  int fake_seq;

};

#endif
/***************************************************************************
 * RCS INFORMATION:
 *
 *	$RCSfile: ComputeNonbondedSelf.h,v $
 *	$Author: jim $	$Locker:  $		$State: Exp $
 *	$Revision: 1.6 $	$Date: 1997/01/16 20:00:21 $
 *
 ***************************************************************************
 * REVISION HISTORY:
 *
 * $Log: ComputeNonbondedSelf.h,v $
 * Revision 1.6  1997/01/16 20:00:21  jim
 * Added reduction calls to ComputeNonbondedSelf and ...Pair.
 * Also moved some code from ...Excl to ...Util.
 *
 * Revision 1.5  1996/11/05 22:40:05  jim
 * commented out undefined virtual destructor
 *
 * Revision 1.4  1996/10/31 21:57:41  jim
 * first incarnation as ComputeNonbondedSelf
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

