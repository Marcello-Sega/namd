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

#include "Sequencer.h"
#include "HomePatch.h"

#define MIN_DEBUG_LEVEL 4
#define DEBUGM
#include "Debug.h"

void SequencerThreadRun(Sequencer* arg)
{
    arg->threadRun();
}

void Sequencer::run(int numberOfCycles)
{
    this->numberOfCycles = numberOfCycles;
    thread = CthCreate((CthVoidFn)&(SequencerThreadRun),(void*)(this),0);
    CthSetStrategyDefault(thread);
    CthAwaken(thread);
}

void Sequencer::threadRun(void)
{
    const int stepsPerCycle = 2;
    const int timestep = 1.0;
    int step, cycle;
    for ( cycle = 0; cycle < numberOfCycles; ++cycle )
    {
        for ( step = 0; step < stepsPerCycle; ++step )
        {
            // patch->addForceToMomentum(0.5*timestep);
            // patch->addVelocityToPosition(timestep);
	    DebugM(3,"(" << cycle << "," << step << ") "
		<< "Sending positionsReady().\n");
            patch->positionsReady();
	    DebugM(2,"(" << cycle << "," << step << ") "
		<< "Suspending.\n");
            suspend();
	    DebugM(2,"(" << cycle << "," << step << ") "
		<< "Awakened!\n");
            // patch->addForceToMomentum(0.5*timestep);
        }
    }
    DebugM(4,"Exiting.\n");
    terminate();
}
