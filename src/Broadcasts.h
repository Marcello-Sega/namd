//-*-c++-*-
/***************************************************************************/
/*                                                                         */
/*              (C) Copyright 1996 The Board of Trustees of the            */
/*                          University of Illinois                         */
/*                           All Rights Reserved                           */
/*									   */
/***************************************************************************/

/***************************************************************************
 * DESCRIPTION:
 *
 ***************************************************************************/

#ifndef BROADCASTS_H
#define BROADCASTS_H

#include "NamdTypes.h"
#include "Lattice.h"
#include "BroadcastObject.h"

// Tags used in common by all users of broadcast system.
enum {
  velocityRescaleFactorTag,	// used in velocity rescaling
  periodicLatticeTag,	// used in constant pressure to send lattice
  dummyTag
};

// Broadcasts used by Contoller <-> Sequencer communication.
struct ControllerBroadcasts
{
  SimpleBroadcastObject<BigReal> velocityRescaleFactor;
  SimpleBroadcastObject<BigReal> positionRescaleFactor;
  // SimpleBroadcastObject<Lattice> lattice;
#ifdef CYCLE_BARRIER
  SimpleBroadcastObject<int> cycleBarrier;
#endif
};

#endif // BROADCASTS_H

/***************************************************************************
 * RCS INFORMATION:
 *
 *	$RCSfile $
 *	$Author $	$Locker:  $		$State: Exp $
 *	$Revision: 1.2 $	$Date: 1997/08/22 19:27:34 $
 *
 ***************************************************************************
 * REVISION HISTORY:
 *
 * $Log: Broadcasts.h,v $
 * Revision 1.2  1997/08/22 19:27:34  brunner
 * Added cycle barrier, enabled by compiling with -DCYCLE_BARRIER
 *
 * Revision 1.1  1997/03/21 23:05:31  jim
 * Added Berendsen's pressure coupling method, won't work with MTS yet.
 *
 *
 ***************************************************************************/
