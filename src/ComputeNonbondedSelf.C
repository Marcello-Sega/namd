/***************************************************************************/
/*                                                                         */
/*              (C) Copyright 1996 The Board of Trustees of the            */
/*                          University of Illinois                         */
/*                           All Rights Reserved                           */
/*									   */
/***************************************************************************/

/***************************************************************************
 * DESCRIPTION: ComputeNonbondedSelf.C
 *
 ***************************************************************************/

#include "ComputeNonbondedSelf.h"
#include "ReductionMgr.h"

#define MIN_DEBUG_LEVEL 5
#define DEBUGM
#include "Debug.h"

ComputeNonbondedSelf::ComputeNonbondedSelf(ComputeID c, PatchID pid)
  : ComputePatch(c,pid)
{
  reduction = ReductionMgr::Object();
  registerReductionData(reduction);
  fake_seq = 0;
}


ComputeNonbondedSelf::~ComputeNonbondedSelf()
{
  unregisterReductionData(reduction);
}


void ComputeNonbondedSelf::doForce(Position* p,
                               Force* f,
                               AtomProperties* a)
{
  DebugM(2,"doForce() called.\n");
  DebugM(1,numAtoms << " patch 1 atoms\n");

  BigReal reductionData[reductionDataSize];
  for ( int i = 0; i < reductionDataSize; ++i ) reductionData[i] = 0;

  ComputeNonbondedUtil::calcSelf(p,f,a,numAtoms,reductionData);

  submitReductionData(reductionData,reduction,fake_seq);
  ++fake_seq;
}

/***************************************************************************
 * RCS INFORMATION:
 *
 *	$RCSfile: ComputeNonbondedSelf.C,v $
 *	$Author: ari $	$Locker:  $		$State: Exp $
 *	$Revision: 1.777 $	$Date: 1997/01/17 19:35:58 $
 *
 ***************************************************************************
 * REVISION HISTORY:
 *
 * $Log: ComputeNonbondedSelf.C,v $
 * Revision 1.777  1997/01/17 19:35:58  ari
 * Internal CVS leveling release.  Start development code work
 * at 1.777.1.1.
 *
 * Revision 1.10  1997/01/16 20:00:16  jim
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
 * Revision 1.6  1996/11/05 21:12:45  jim
 * fixed modified pairs
 *
 * Revision 1.5  1996/11/05 05:08:56  jim
 * added nonbonded compute code for one case ( no ifs )
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

