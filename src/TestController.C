/**
***  Copyright (c) 1995, 1996, 1997, 1998, 1999, 2000 by
***  The Board of Trustees of the University of Illinois.
***  All rights reserved.
**/

#include "Node.h"
#include "Molecule.h"
#include "SimParameters.h"
#include "TestController.h"
#include "ReductionMgr.h"
#include "CollectionMaster.h"
#include "Output.h"
#include "strlib.h"
#include "BroadcastObject.h"
#include "NamdState.h"
#include "Broadcasts.h"
#include "LdbCoordinator.h"
#include "Thread.h"
#include <math.h>

#ifndef cbrt
  // cbrt() not in math.h on goneril
  #define cbrt(x)  pow(x,(double)(1.0/3.0))
#endif

#define MIN_DEBUG_LEVEL 3
//#define DEBUGM
#include "Debug.h"

TestController::TestController(NamdState *s) :
	Controller(s)
{
  ;
}

TestController::~TestController(void)
{
  ;
}

extern int eventEndOfTimeStep;

void TestController::algorithm(int)
{
    int step = simParams->firstTimestep;

    const int numberOfSteps = simParams->N;
    // const int stepsPerCycle = simParams->stepsPerCycle;
    // const BigReal timestep = simParams->dt;

    iout << iINFO << "**********************************\n";
    iout << iINFO << "* NAMD RUNNING IN SELF-TEST MODE *\n";
    iout << iINFO << "**********************************\n" << endi;

    for ( ; step <= numberOfSteps; ++step )
    {
        enqueueCollections(step);
        traceUserEvent(eventEndOfTimeStep);
        receivePressure(step);
        printEnergies(step);
        rescaleVelocities(step);
	tcoupleVelocities(step);
	berendsenPressure(step);
#ifdef CYCLE_BARRIER
	if (!((step+1) % stepsPerCycle))
	{
	  broadcast->cycleBarrier.publish(step,1);
	  CkPrintf("Cycle time at sync Wall: %f CPU %f\n",
		  CmiWallTimer(),CmiTimer());
	}
#endif
/*
	if ( LdbCoordinator::Object()->balanceNow(step) ) {
	  LdbCoordinator::Object()->rebalance(this);
	}
*/

    }

    terminate();
}

void TestController::berendsenPressure(int step)
{
  const int freq = simParams->berendsenPressureFreq;
  if ( simParams->berendsenPressureOn && !(step%freq) )
  {
    BigReal factor = 0.0001;
    /*
    BigReal factor = pressure - simParams->berendsenPressureTarget;
    factor *= simParams->berendsenPressureCompressibility;
    factor *= ( simParams->dt * freq );
    factor /= simParams->berendsenPressureRelaxationTime;
    */
    factor += 1.0;
    factor = cbrt(factor);
    broadcast->positionRescaleFactor.publish(step,Tensor::identity()*factor);
    state->lattice.rescale(Tensor::identity()*factor);
  }
}

