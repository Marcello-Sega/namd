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
extern "C" void trace_user_event(int event);

void TestController::algorithm(void)
{
    int step = simParams->firstTimestep;

    const int numberOfSteps = simParams->N;
    const int stepsPerCycle = simParams->stepsPerCycle;
    const BigReal timestep = simParams->dt;

    iout << iINFO << "**********************************\n";
    iout << iINFO << "* NAMD RUNNING IN SELF-TEST MODE *\n";
    iout << iINFO << "**********************************\n" << endi;

    for ( ; step <= numberOfSteps; ++step )
    {
        enqueueCollections(step);
        trace_user_event(eventEndOfTimeStep);
        receivePressure(step);
        printEnergies(step);
        rescaleVelocities(step);
	tcoupleVelocities(step);
	berendsenPressure(step);
#ifdef CYCLE_BARRIER
	if (!((step+1) % stepsPerCycle))
	{
	  broadcast->cycleBarrier.publish(step,1);
	  CPrintf("Cycle time at sync Wall: %f CPU %f\n",
		  CmiWallTimer(),CmiTimer());
	}
#endif
	if ( LdbCoordinator::Object()->balanceNow(step) ) {
	  LdbCoordinator::Object()->rebalance(this);
	}

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
    broadcast->positionRescaleFactor.publish(step,factor);
    state->lattice.rescale(factor);
  }
}

/***************************************************************************
 * RCS INFORMATION:
 *
 *	$RCSfile $
 *	$Author $	$Locker:  $		$State: Exp $
 *	$Revision: 1.4 $	$Date: 1998/09/15 03:06:00 $
 *
 ***************************************************************************
 * REVISION HISTORY:
 *
 * $Log: TestController.C,v $
 * Revision 1.4  1998/09/15 03:06:00  jim
 * Fixed test mode.
 *
 * Revision 1.3  1998/08/11 16:30:31  jim
 * Modified output from periodic boundary simulations to return atoms to
 * internally consistent coordinates.  We store the transformations which
 * were performed and undo them at the end.  It might be better to do this
 * by always keeping the original coordinates and only doing the transform
 * for the nonbonded terms but this works for now.
 *
 * Revision 1.2  1998/04/06 16:34:12  jim
 * Added DPME (single processor only), test mode, and momenta printing.
 *
 * Revision 1.1  1998/03/31 04:55:48  jim
 * Added test mode, fixed errors in virial with full electrostatics.
 *
 *
 ***************************************************************************/
