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
#include "LdbCoordinator.h"

#define MIN_DEBUG_LEVEL 4
// #define DEBUGM
#include "Debug.h"

ComputeNonbondedSelf::ComputeNonbondedSelf(ComputeID c, PatchID pid,
		int minPartition, int maxPartition, int numPartitions)
  : ComputePatch(c,pid),
    minPart(minPartition), maxPart(maxPartition), numParts(numPartitions)
{
  reduction = ReductionMgr::Object()->willSubmit(REDUCTIONS_BASIC);
}

void ComputeNonbondedSelf::initialize() {
  ComputePatch::initialize();
  avgPositionBox = patch->registerAvgPositionPickup(cid);
}

ComputeNonbondedSelf::~ComputeNonbondedSelf()
{
  delete reduction;
  if (avgPositionBox != NULL) {
    patch->unregisterAvgPositionPickup(cid,&avgPositionBox);
  }
}


void ComputeNonbondedSelf::doForce(Position* p,
                               Results* r,
                               AtomProperties* a)
{
  // Inform load balancer. 
  // I assume no threads will suspend until endWork is called
  LdbCoordinator::Object()->startWork(cid,0); // Timestep not used

  DebugM(2,"doForce() called.\n");
  DebugM(1,numAtoms << " patch 1 atoms\n");
  DebugM(3, "NUMATOMSxNUMATOMS = " << numAtoms*numAtoms << "\n");

  BigReal reductionData[reductionDataSize];
  for ( int i = 0; i < reductionDataSize; ++i ) reductionData[i] = 0;

  if ( patch->flags.doNonbonded )
  {
    nonbonded params;
    params.p[0] = p;
    params.p[1] = p;
    params.ff[0] = r->f[Results::nbond];
    params.ff[1] = r->f[Results::nbond];
    params.a[0] = a;
    params.a[1] = a;
    params.numAtoms[0] = numAtoms-1;
    params.numAtoms[1] = numAtoms;
    params.reduction = reductionData;

    params.minPart = minPart;
    params.maxPart = maxPart;
    params.numParts = numParts;

    if ( patch->flags.doFullElectrostatics )
    {
      params.fullf[0] = r->f[Results::slow];
      params.fullf[1] = r->f[Results::slow];
      if ( patch->flags.doMolly ) {
        calcSelf(&params);
        Position *p_avg = avgPositionBox->open();
        params.p[0] = p_avg;
        params.p[1] = p_avg;
        calcSlowSelf(&params);
        avgPositionBox->close(&p_avg);
      } else {
        calcFullSelf(&params);
      }
    }
    else
      calcSelf(&params);
  }

  submitReductionData(reductionData,reduction);
  reduction->submit();
  // Inform load balancer
  LdbCoordinator::Object()->endWork(cid,0); // Timestep not used
}

/***************************************************************************
 * RCS INFORMATION:
 *
 *	$RCSfile: ComputeNonbondedSelf.C,v $
 *	$Author: jim $	$Locker:  $		$State: Exp $
 *	$Revision: 1.1015 $	$Date: 1999/08/20 19:11:10 $
 *
 ***************************************************************************
 * REVISION HISTORY:
 *
 * $Log: ComputeNonbondedSelf.C,v $
 * Revision 1.1015  1999/08/20 19:11:10  jim
 * Added MOLLY - mollified impluse method.
 *
 * Revision 1.1014  1999/06/17 17:05:40  jim
 * Renamed seq to step in most places.  Now has meaning only to user.
 *
 * Revision 1.1013  1999/06/17 15:46:10  jim
 * Completely rewrote reduction system to eliminate need for sequence numbers.
 *
 * Revision 1.1012  1999/03/18 02:41:16  jim
 * Turned off stray DEBUGM code.
 *
 * Revision 1.1011  1998/07/02 21:06:36  jim
 * Added support for splitting ComputeNonbondedSelf into multiple computes.
 *
 * Revision 1.1010  1997/05/20 15:49:10  nealk
 * Pair, Self, and Excl not use the same parameters!
 *
 * Revision 1.1009  1997/05/15 17:43:48  nealk
 * Merged Pair and Self to use same headers.
 *
 * Revision 1.1008  1997/04/08 07:08:25  ari
 * Modification for dynamic loadbalancing - moving computes
 * Still bug in new computes or usage of proxies/homepatches.
 * Works if ldbStrategy is none as before.
 *
 * Revision 1.1007  1997/04/06 22:45:02  ari
 * Add priorities to messages.  Mods to help proxies without computes.
 * Added quick enhancement to end of list insertion of ResizeArray(s)
 *
 * Revision 1.1006  1997/03/27 20:25:43  brunner
 * Changes for LdbCoordinator, the load balance control BOC
 *
 * Revision 1.1005  1997/03/25 23:00:59  jim
 * Added nonbondedFrequency parameter and multiple time-stepping
 *
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

