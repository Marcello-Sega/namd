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
};

#endif // BROADCASTS_H

/***************************************************************************
 * RCS INFORMATION:
 *
 *	$RCSfile $
 *	$Author $	$Locker:  $		$State: Exp $
 *	$Revision: 1.1 $	$Date: 1997/03/21 23:05:31 $
 *
 ***************************************************************************
 * REVISION HISTORY:
 *
 * $Log: Broadcasts.h,v $
 * Revision 1.1  1997/03/21 23:05:31  jim
 * Added Berendsen's pressure coupling method, won't work with MTS yet.
 *
 *
 ***************************************************************************/
