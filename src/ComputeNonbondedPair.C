/***************************************************************************/
/*                                                                         */
/*              (C) Copyright 1996 The Board of Trustees of the            */
/*                          University of Illinois                         */
/*                           All Rights Reserved                           */
/*									   */
/***************************************************************************/

/***************************************************************************
 * DESCRIPTION: ComputeNonbondedPair.C
 *
 ***************************************************************************/

#include "ComputeNonbondedPair.h"
#include "ReductionMgr.h"

#define MIN_DEBUG_LEVEL 5
#define DEBUGM
#include "Debug.h"

ComputeNonbondedPair::ComputeNonbondedPair(ComputeID c, PatchID pid[])
  : ComputePatchPair(c,pid)
{
  reduction = ReductionMgr::Object();
  registerReductionData(reduction);
  fake_seq = 0;
}


ComputeNonbondedPair::~ComputeNonbondedPair()
{
  unregisterReductionData(reduction);
}


void ComputeNonbondedPair::doForce(Position* p[2],
                               Force* f[2],
                               AtomProperties* a[2])
{
  DebugM(2,"doForce() called.\n");
  DebugM(1, numAtoms[0] << " patch 1 atoms and " <<
	numAtoms[1] << " patch 2 atoms\n");

  BigReal reductionData[reductionDataSize];
  for ( int i = 0; i < reductionDataSize; ++i ) reductionData[i] = 0;

  ComputeNonbondedUtil::calcPair(p,f,a,numAtoms,reductionData);

  submitReductionData(reductionData,reduction,fake_seq);
  ++fake_seq;
}

/***************************************************************************
 * RCS INFORMATION:
 *
 *	$RCSfile: ComputeNonbondedPair.C,v $
 *	$Author: jim $	$Locker:  $		$State: Exp $
 *	$Revision: 1.10 $	$Date: 1997/01/16 20:00:06 $
 *
 ***************************************************************************
 * REVISION HISTORY:
 *
 * $Log: ComputeNonbondedPair.C,v $
 * Revision 1.10  1997/01/16 20:00:06  jim
 * Added reduction calls to ComputeNonbondedSelf and ...Pair.
 * Also moved some code from ...Excl to ...Util.
 *
 * Revision 1.9  1996/11/30 21:09:15  jim
 * cleaned up debug messages
 *
 * Revision 1.8  1996/11/21 00:59:16  jim
 * moved ComputeNonbondedUtil::select() call to SimParameters
 *
 * Revision 1.7  1996/11/21 00:22:01  jim
 * first time using ComputeNonbondedUtil functions
 *
 * Revision 1.6  1996/11/05 21:12:12  jim
 * fixed modified pairs
 *
 * Revision 1.5  1996/11/05 05:08:56  jim
 * added nonbonded compute code for one case ( no ifs )
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

