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

#define MIN_DEBUG_LEVEL 5
//#define DEBUGM
#include "Debug.h"

ComputeNonbondedPair::ComputeNonbondedPair(ComputeID c, PatchID pid[], int trans[])
  : ComputePatchPair(c,pid,trans)
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
  DebugM(5, numAtoms[0] << " patch #1 atoms and " <<
	numAtoms[1] << " patch #2 atoms\n");

  BigReal reductionData[reductionDataSize];
  for ( int i = 0; i < reductionDataSize; ++i ) reductionData[i] = 0;

  if ( numAtoms[0] && numAtoms[1] )
  {
    // swap to place more atoms in inner loop (second patch)
    if ( numAtoms[0] > numAtoms[1] )
    {
      Position* p_r[2];
      p_r[0] = p[1];
      p_r[1] = p[0];
      Force* f_r[2];
      f_r[0] = f[1];
      f_r[1] = f[0];
      AtomProperties* a_r[2];
      a_r[0] = a[1];
      a_r[1] = a[0];
      int numAtoms_r[2];
      numAtoms_r[0] = numAtoms[1];
      numAtoms_r[1] = numAtoms[0];
      if ( patch[0]->flags.doFullElectrostatics )
        ComputeNonbondedUtil::calcFullPair(p_r,f_r,f_r,a_r,numAtoms_r,reductionData);
      else
        ComputeNonbondedUtil::calcPair(p_r,f_r,a_r,numAtoms_r,reductionData);
    }
    else
    {
      if ( patch[0]->flags.doFullElectrostatics )
        ComputeNonbondedUtil::calcFullPair(p,f,f,a,numAtoms,reductionData);
      else
        ComputeNonbondedUtil::calcPair(p,f,a,numAtoms,reductionData);
    }
  }

  submitReductionData(reductionData,reduction,fake_seq);
  ++fake_seq;
}

/***************************************************************************
 * RCS INFORMATION:
 *
 *	$RCSfile: ComputeNonbondedPair.C,v $
 *	$Author: jim $	$Locker:  $		$State: Exp $
 *	$Revision: 1.1002 $	$Date: 1997/02/28 04:47:04 $
 *
 ***************************************************************************
 * REVISION HISTORY:
 *
 * $Log: ComputeNonbondedPair.C,v $
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

