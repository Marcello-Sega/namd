/**
***  Copyright (c) 1995, 1996, 1997, 1998, 1999, 2000 by
***  The Board of Trustees of the University of Illinois.
***  All rights reserved.
**/

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

