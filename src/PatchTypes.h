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

#ifndef PATCHTYPES_H
#define PATCHTYPES_H

#include "NamdTypes.h"
#include "Lattice.h"

class Flags
{
public:
  int step;			// timestep number reported to user
  				// Same number may appear multiple times!
  int doNonbonded;
  int doFullElectrostatics;
  int submitLoadStats;
  int maxForceUsed;		// may ignore slower force classes
  int maxForceMerged;		// add this and faster to normal

  Lattice lattice;		// rather than shipping around separately

};

class Results
{
public:
  enum { normal=0, nbond=1, slow=2, maxNumForces=3 };
  Force *f[maxNumForces];
};

#endif


/***************************************************************************
 * RCS INFORMATION:
 *
 *	$RCSfile $
 *	$Author $	$Locker:  $		$State: Exp $
 *	$Revision: 1.10 $	$Date: 1999/06/17 17:05:46 $
 *
 ***************************************************************************
 * REVISION HISTORY:
 *
 * $Log: PatchTypes.h,v $
 * Revision 1.10  1999/06/17 17:05:46  jim
 * Renamed seq to step in most places.  Now has meaning only to user.
 *
 * Revision 1.9  1997/12/22 21:29:25  jim
 * Proxies no longer send empty arrays back to HomePatch.  Requires some new
 * flags to be set correctly in Sequencer in order to work.  These are:
 *   maxForceMerged - this and faster are added into Results::normal array
 *   maxForceUsed - all forces slower than this are discarded (assumed zero)
 * Generally maxForceMerged doesn't change but maxForceUsed depends on timestep.
 *
 * Revision 1.8  1997/04/10 09:14:05  ari
 * Final debugging for compute migration / proxy creation for load balancing.
 * Lots of debug code added, mostly turned off now.
 * Fixed bug in PositionBox when Patch had no dependencies.
 * Eliminated use of cout and misuse of iout in numerous places.
 *                                            Ari & Jim
 *
 * Revision 1.7  1997/04/08 21:08:45  jim
 * Contant pressure now correct on multiple nodes, should work with MTS.
 *
 * Revision 1.6  1997/03/27 20:25:50  brunner
 * Changes for LdbCoordinator, the load balance control BOC
 *
 * Revision 1.5  1997/03/25 23:01:00  jim
 * Added nonbondedFrequency parameter and multiple time-stepping
 *
 * Revision 1.4  1997/03/19 11:54:48  ari
 * Add Broadcast mechanism.
 * Fixed RCS Log entries on files that did not have Log entries.
 * Added some register variables to Molecule and ComputeNonbondedExcl.C
 *
 ***************************************************************************/
