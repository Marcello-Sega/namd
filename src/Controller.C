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
#include "ReductionMgr.h"

#define MIN_DEBUG_LEVEL 4
#define DEBUGM
#include "Debug.h"

Controller::Controller(NamdState *s) :
	state(s),
	simParams(Node::Object()->simParameters),
	reduction(ReductionMgr::Object())
{
    reduction->subscribe(REDUCTION_KINETIC_ENERGY);
}

Controller::~Controller(void)
{
    reduction->unsubscribe(REDUCTION_KINETIC_ENERGY);
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
    int seq = 0;
    BigReal ke;
    reduction->require(seq, REDUCTION_KINETIC_ENERGY, ke);
    DebugM(5, "Step: " << seq << " KE: " << ke << "\n");
    ++seq;
    for ( cycle = 0; cycle < numberOfCycles; ++cycle )
    {
        for ( step = 0; step < stepsPerCycle; ++step )
        {
	    reduction->require(seq, REDUCTION_KINETIC_ENERGY, ke);
	    DebugM(5, "Step: " << seq << " KE: " << ke << "\n");
	    ++seq;
        }
    }
    DebugM(4, "Controller: Exiting.\n");
    terminate();
}
