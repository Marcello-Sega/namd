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
#include "SimParameters.h"
#include "Sequencer.h"
#include "HomePatch.h"

#define MIN_DEBUG_LEVEL 4
#define DEBUGM
#include "Debug.h"

Sequencer::Sequencer(HomePatch *p) :
	patch(p),
	simParams(Node::Object()->simParameters)
{
  ;
}

void Sequencer::threadRun(Sequencer* arg)
{
    arg->algorithm();
}

void Sequencer::run(int numberOfCycles)
{
    stepsPerCycle = simParams->stepsPerCycle;
    if ( numberOfCycles ) this->numberOfCycles = numberOfCycles;
    else this->numberOfCycles = simParams->N / stepsPerCycle;
    thread = CthCreate((CthVoidFn)&(threadRun),(void*)(this),0);
    CthSetStrategyDefault(thread);
    CthAwaken(thread);
}

void Sequencer::algorithm(void)
{
    const int numberOfCycles = this->numberOfCycles;
    const int stepsPerCycle = this->stepsPerCycle;
    const BigReal timestep = simParams->dt;
    int step, cycle;
    patch->positionsReady();
    suspend();
    for ( cycle = 0; cycle < numberOfCycles; ++cycle )
    {
        for ( step = 0; step < stepsPerCycle; ++step )
        {
            patch->addForceToMomentum(0.5*timestep);
            patch->addVelocityToPosition(timestep);
	    DebugM(3, patch->getPatchID()
		<< ": (" << cycle << "," << step << ") "
		<< "Sending positionsReady().\n");
            patch->positionsReady();
	    DebugM(2, patch->getPatchID()
		<< ": (" << cycle << "," << step << ") "
		<< "Suspending.\n");
            suspend();
	    DebugM(2, patch->getPatchID()
		<< ": (" << cycle << "," << step << ") "
		<< "Awakened!\n");
            patch->addForceToMomentum(0.5*timestep);
        }
    }
    DebugM(4, patch->getPatchID() << ": Exiting.\n");
    terminate();
}
