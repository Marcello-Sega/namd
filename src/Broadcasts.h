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
  SimpleBroadcastObject<Vector> positionRescaleFactor;
  SimpleBroadcastObject<BigReal> tcoupleCoefficient;
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
 *	$Revision: 1.4 $	$Date: 1999/01/06 19:19:19 $
 *
 ***************************************************************************
 * REVISION HISTORY:
 *
 * $Log: Broadcasts.h,v $
 * Revision 1.4  1999/01/06 19:19:19  jim
 * Broadcast and Sequencers understand anisotropic volume rescaling factors.
 *
 * Revision 1.3  1998/03/06 20:55:24  jim
 * Added temperature coupling.
 *
 * Revision 1.2  1997/08/22 19:27:34  brunner
 * Added cycle barrier, enabled by compiling with -DCYCLE_BARRIER
 *
 * Revision 1.1  1997/03/21 23:05:31  jim
 * Added Berendsen's pressure coupling method, won't work with MTS yet.
 *
 *
 ***************************************************************************/
