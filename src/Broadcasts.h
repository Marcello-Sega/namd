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

#define SCRIPT_END 0
#define SCRIPT_RUN 1
#define SCRIPT_OUTPUT 2
#define SCRIPT_MEASURE 3
#define SCRIPT_REINITVELS 4
#define SCRIPT_CHECKPOINT 5
#define SCRIPT_REVERT 6
#define SCRIPT_MINIMIZE 7

// Tags used in common by all users of broadcast system.
enum {
  velocityRescaleFactorTag,
  positionRescaleFactorTag,
  tcoupleCoefficientTag,
  minimizeCoefficientTag,
#if USE_BARRIER
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
  SimpleBroadcastObject<BigReal> minimizeCoefficient;
#if USE_BARRIER
  SimpleBroadcastObject<int> cycleBarrier;
#endif
  SimpleBroadcastObject<int> scriptBarrier;

  ControllerBroadcasts() : 
    velocityRescaleFactor(velocityRescaleFactorTag),
    positionRescaleFactor(positionRescaleFactorTag),
    tcoupleCoefficient(tcoupleCoefficientTag),
    minimizeCoefficient(minimizeCoefficientTag),
#if USE_BARRIER
    cycleBarrier(cycleBarrierTag),
#endif
    scriptBarrier(scriptBarrierTag)
  { ; }
};

#endif // BROADCASTS_H

