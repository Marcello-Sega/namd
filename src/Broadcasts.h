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
  velocityRescaleFactorTag,
  positionRescaleFactorTag,
  tcoupleCoefficientTag,
#ifdef CYCLE_BARRIER
  cycleBarrierTag,
#endif
  scriptBarrierTag,
  dummyTag
};

// Broadcasts used by Contoller <-> Sequencer communication.
struct ControllerBroadcasts
{
  SimpleBroadcastObject<BigReal> velocityRescaleFactor;
  SimpleBroadcastObject<Tensor> positionRescaleFactor;
  SimpleBroadcastObject<BigReal> tcoupleCoefficient;
#ifdef CYCLE_BARRIER
  SimpleBroadcastObject<int> cycleBarrier;
#endif
  SimpleBroadcastObject<int> scriptBarrier;

  ControllerBroadcasts() : 
    velocityRescaleFactor(velocityRescaleFactorTag),
    positionRescaleFactor(positionRescaleFactorTag),
    tcoupleCoefficient(tcoupleCoefficientTag),
#ifdef CYCLE_BARRIER
    cycleBarrier(cycleBarrierTag),
#endif
    scriptBarrier(scriptBarrierTag)
  { ; }
};

#endif // BROADCASTS_H

/***************************************************************************
 * RCS INFORMATION:
 *
 *	$RCSfile $
 *	$Author $	$Locker:  $		$State: Exp $
 *	$Revision: 1.6 $	$Date: 1999/09/03 20:46:05 $
 *
 ***************************************************************************
 * REVISION HISTORY:
 *
 * $Log: Broadcasts.h,v $
 * Revision 1.6  1999/09/03 20:46:05  jim
 * Support for non-orthogonal periodic boundary conditions.
 *
 * Revision 1.5  1999/05/26 22:23:53  jim
 * Added basic Tcl scripting, fixed bugs in broadcasts.
 *
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
