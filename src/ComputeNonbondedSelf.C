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
#include "Patch.h"

#define MIN_DEBUG_LEVEL 4
#define DEBUGM
#include "Debug.h"

ComputeNonbondedSelf::ComputeNonbondedSelf(ComputeID c, PatchID pid)
  : ComputePatch(c,pid)
{
  reduction = ReductionMgr::Object();
  registerReductionData(reduction);
}


ComputeNonbondedSelf::~ComputeNonbondedSelf()
{
  unregisterReductionData(reduction);
}


void ComputeNonbondedSelf::doForce(Position* p,
                               Results* r,
                               AtomProperties* a)
{
  DebugM(2,"doForce() called.\n");
  DebugM(1,numAtoms << " patch 1 atoms\n");
  DebugM(3, "NUMATOMSxNUMATOMS = " << numAtoms*numAtoms << "\n");

  BigReal reductionData[reductionDataSize];
  for ( int i = 0; i < reductionDataSize; ++i ) reductionData[i] = 0;

  if ( patch->flags.doFullElectrostatics )
    calcFullSelf(p,r->f[Results::normal],r->f[Results::slow],
	a,numAtoms,reductionData);
  else
    calcSelf(p,r->f[Results::normal],a,numAtoms,reductionData);

  submitReductionData(reductionData,reduction,patch->flags.seq);
}

/***************************************************************************
 * RCS INFORMATION:
 *
 *	$RCSfile: ComputeNonbondedSelf.C,v $
 *	$Author: jim $	$Locker:  $		$State: Exp $
 *	$Revision: 1.1004 $	$Date: 1997/03/18 21:35:30 $
 *
 ***************************************************************************
 * REVISION HISTORY:
 *
 * $Log: ComputeNonbondedSelf.C,v $
 * Revision 1.1004  1997/03/18 21:35:30  jim
 * Eliminated fake_seq.  Reductions now use Patch::flags.seq.
 *
 * Revision 1.1003  1997/03/13 06:37:03  jim
 * Multiple time-stepping implemented, still needs proper splitting functions.
 *
 * Revision 1.1002  1997/03/10 17:40:10  ari
 * UniqueSet changes - some more commenting and cleanup
 *
 * Revision 1.1001  1997/02/28 04:47:05  jim
 * Full electrostatics now works with fulldirect on one node.
 *
 * Revision 1.1000  1997/02/06 15:58:12  ari
 * Resetting CVS to merge branches back into the main trunk.
 * We will stick to main trunk development as suggested by CVS manual.
 * We will set up tags to track fixed points of development/release
 * as suggested by CVS manual - all praise the CVS manual.
 *
 * Revision 1.778  1997/01/28 00:30:24  ari
 * internal release uplevel to 1.778
 *
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

