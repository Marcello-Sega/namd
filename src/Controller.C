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
#include "Controller.h"

#define MIN_DEBUG_LEVEL 4
#define DEBUGM
#include "Debug.h"

Controller::Controller(NamdState *s) :
	state(s),
	simParams(Node::Object()->simParameters),
	reduction(ReductionMgr::Object())
{
  ;
}

void Controller::threadRun(Controller* arg)
{
    arg->algorithm();
}

void Controller::run(int numberOfCycles)
{
    stepsPerCycle = simParams->stepsPerCycle;
    if ( numberOfCycles ) this->numberOfCycles = numberOfCycles;
    else this->numberOfCycles = simParams->N / stepsPerCycle;
    thread = CthCreate((CthVoidFn)&(threadRun),(void*)(this),0);
    CthSetStrategyDefault(thread);
    CthAwaken(thread);
}

void Controller::algorithm(void)
{
    DebugM(4, "Controller algorithm active.\n");
    const int numberOfCycles = this->numberOfCycles;
    const int stepsPerCycle = this->stepsPerCycle;
    const BigReal timestep = simParams->dt;
    int step, cycle;
    for ( cycle = 0; cycle < numberOfCycles; ++cycle )
    {
        for ( step = 0; step < stepsPerCycle; ++step )
        {
	    DebugM(2, patch->getPatchID()
		<< ": (" << cycle << "," << step << ") "
		<< "Awakened!\n");
        }
    }
    DebugM(4, "Controller: Exiting.\n");
    terminate();
}
