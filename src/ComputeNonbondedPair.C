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
#include "Patch.h"
#include "LdbCoordinator.h"
#include "Priorities.h"
#include "PatchMap.h"

#define MIN_DEBUG_LEVEL 4
#define DEBUGM
#include "Debug.h"

ComputeNonbondedPair::ComputeNonbondedPair(ComputeID c, PatchID pid[], int trans[])
  : ComputePatchPair(c,pid,trans)
{
  reduction = ReductionMgr::Object();
  registerReductionData(reduction);
}


ComputeNonbondedPair::~ComputeNonbondedPair()
{
  unregisterReductionData(reduction);
}


int ComputeNonbondedPair::priority(void)
{
  PatchMap *patchMap = PatchMap::Object();
  int myNode = CMyPe();
  if ( patchMap->node(patchID[0]) != myNode ||
       patchMap->node(patchID[1]) != myNode )
  {
    return Priorities::nonlocal;
  }
  else
  {
    return Priorities::base;
  }
}


int ComputeNonbondedPair::noWork() {

  // return 0;  // for testing
  if ( numAtoms[0] && numAtoms[1] && patch[0]->flags.doNonbonded )
  {
    return 0;  // work to do, enqueue as usual
  } else
  {
    // Inform load balancer
    LdbCoordinator::Object()->startWork(cid,0); // Timestep not used
    // fake out patches and reduction system

    BigReal reductionData[reductionDataSize];
    for ( int i = 0; i < reductionDataSize; ++i ) reductionData[i] = 0;

    Position* p[2];
    Results* r[2];
    AtomProperties* a[2];

    // Open up positionBox, forceBox, and atomBox
    for (i=0; i<2; i++) {
      p[i] = positionBox[i]->open();
      r[i] = forceBox[i]->open();
      a[i] = atomBox[i]->open();
    }

    // Close up boxes
    for (i=0; i<2; i++) {
      positionBox[i]->close(&p[i]);
      forceBox[i]->close(&r[i]);
      atomBox[i]->close(&a[i]);
    }

    submitReductionData(reductionData,reduction,patch[0]->flags.seq);

    // Inform load balancer
    LdbCoordinator::Object()->endWork(cid,0); // Timestep not used

    return 1;  // no work to do, do not enqueue
  }
}


void ComputeNonbondedPair::doForce(Position* p[2],
                               Results* r[2],
                               AtomProperties* a[2])
{
  // Inform load balancer. 
  // I assume no threads will suspend until endWork is called
  LdbCoordinator::Object()->startWork(cid,0); // Timestep not used

  DebugM(2,"doForce() called.\n");
  DebugM(2, numAtoms[0] << " patch #1 atoms and " <<
	numAtoms[1] << " patch #2 atoms\n");

  BigReal reductionData[reductionDataSize];
  for ( int i = 0; i < reductionDataSize; ++i ) reductionData[i] = 0;

  if ( numAtoms[0] && numAtoms[1] )
  {
    nonbonded params;
    params.reduction = reductionData;

    // swap to place more atoms in inner loop (second patch)
    if ( numAtoms[0] > numAtoms[1] )
    {
      params.p[0] = p[1];
      params.p[1] = p[0];
      params.ff[0] = r[1]->f[Results::nbond];
      params.ff[1] = r[0]->f[Results::nbond];
      params.a[0] = a[1];
      params.a[1] = a[0];
      params.numAtoms[0] = numAtoms[1];
      params.numAtoms[1] = numAtoms[0];
      DebugM(3, "NUMATOMSxNUMATOMS = " << numAtoms_r[0]*numAtoms_r[1] << "\n" );
      if ( patch[0]->flags.doFullElectrostatics )
	{
        params.fullf[0] = r[1]->f[Results::slow];
        params.fullf[1] = r[0]->f[Results::slow];
        calcFullPair(&params);
	}
      else
        calcPair(&params);
    }
    else
    {
      params.p[0] = p[0];
      params.p[1] = p[1];
      params.a[0] = a[0];
      params.a[1] = a[1];
      params.numAtoms[0] = numAtoms[0];
      params.numAtoms[1] = numAtoms[1];
      params.ff[0] = r[0]->f[Results::nbond];
      params.ff[1] = r[1]->f[Results::nbond];
      if ( patch[0]->flags.doFullElectrostatics )
	{
        params.fullf[0] = r[0]->f[Results::slow];
        params.fullf[1] = r[1]->f[Results::slow];
        calcFullPair(&params);
	}
      else
        calcPair(&params);
    }
  }

  submitReductionData(reductionData,reduction,patch[0]->flags.seq);

  // Inform load balancer
  LdbCoordinator::Object()->endWork(cid,0); // Timestep not used
}

/***************************************************************************
 * RCS INFORMATION:
 *
 *	$RCSfile: ComputeNonbondedPair.C,v $
 *	$Author: nealk $	$Locker:  $		$State: Exp $
 *	$Revision: 1.1013 $	$Date: 1997/05/20 15:49:09 $
 *
 ***************************************************************************
 * REVISION HISTORY:
 *
 * $Log: ComputeNonbondedPair.C,v $
 * Revision 1.1013  1997/05/20 15:49:09  nealk
 * Pair, Self, and Excl not use the same parameters!
 *
 * Revision 1.1012  1997/04/10 09:13:53  ari
 * Final debugging for compute migration / proxy creation for load balancing.
 * Lots of debug code added, mostly turned off now.
 * Fixed bug in PositionBox when Patch had no dependencies.
 * Eliminated use of cout and misuse of iout in numerous places.
 *                                            Ari & Jim
 *
 * Revision 1.1011  1997/04/08 07:08:21  ari
 * Modification for dynamic loadbalancing - moving computes
 * Still bug in new computes or usage of proxies/homepatches.
 * Works if ldbStrategy is none as before.
 *
 * Revision 1.1010  1997/04/06 22:45:00  ari
 * Add priorities to messages.  Mods to help proxies without computes.
 * Added quick enhancement to end of list insertion of ResizeArray(s)
 *
 * Revision 1.1009  1997/04/03 23:22:19  jim
 * Added basic priority() method to Compute.  Only distinguishes between
 * local and nonlocal computations for now.
 *
 * Revision 1.1008  1997/03/27 20:25:42  brunner
 * Changes for LdbCoordinator, the load balance control BOC
 *
 * Revision 1.1007  1997/03/25 23:00:57  jim
 * Added nonbondedFrequency parameter and multiple time-stepping
 *
 * Revision 1.1006  1997/03/18 21:35:28  jim
 * Eliminated fake_seq.  Reductions now use Patch::flags.seq.
 *
 * Revision 1.1005  1997/03/13 06:36:59  jim
 * Multiple time-stepping implemented, still needs proper splitting functions.
 *
 * Revision 1.1004  1997/03/12 23:59:41  jim
 * Added Compute::noWork() protocol to not enqueue do-nothing compute objects.
 *
 * Revision 1.1003  1997/03/10 17:40:09  ari
 * UniqueSet changes - some more commenting and cleanup
 *
 * Revision 1.1002  1997/02/28 04:47:04  jim
 * Full electrostatics now works with fulldirect on one node.
 *
 * Revision 1.1001  1997/02/13 23:17:16  ari
 * Fixed a final bug in AtomMigration - numatoms in ComputePatchPair.C not
 * set correctly in atomUpdate()
 *
 * Revision 1.1000  1997/02/06 15:58:10  ari
 * Resetting CVS to merge branches back into the main trunk.
 * We will stick to main trunk development as suggested by CVS manual.
 * We will set up tags to track fixed points of development/release
 * as suggested by CVS manual - all praise the CVS manual.
 *
 * Revision 1.779  1997/02/06 15:53:05  ari
 * Updating Revision Line, getting rid of branches
 *
 * Revision 1.778.2.1  1997/02/06 02:35:19  jim
 * Implemented periodic boundary conditions - may not work with
 * atom migration yet, but doesn't seem to alter calculation,
 * appears to work correctly when turned on.
 * NamdState chdir's to same directory as config file in argument.
 *
 * Revision 1.778  1997/01/28 00:30:22  ari
 * internal release uplevel to 1.778
 *
 * Revision 1.777.2.1  1997/01/23 06:17:34  jim
 * Check for empty patches and swap so larger patch is second.
 *
 * Revision 1.777  1997/01/17 19:35:56  ari
 * Internal CVS leveling release.  Start development code work
 * at 1.777.1.1.
 *
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

