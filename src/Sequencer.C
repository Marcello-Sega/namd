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
#include "ReductionMgr.h"
#include "CollectionMgr.h"

#define MIN_DEBUG_LEVEL 5
#define DEBUGM
#include "Debug.h"

Sequencer::Sequencer(HomePatch *p) :
	patch(p),
	simParams(Node::Object()->simParameters),
	reduction(ReductionMgr::Object()),
	collection(CollectionMgr::Object())
{
    reduction->Register(REDUCTION_KINETIC_ENERGY);
}

Sequencer::~Sequencer(void)
{
    reduction->unRegister(REDUCTION_KINETIC_ENERGY);
}

void Sequencer::threadRun(Sequencer* arg)
{
    arg->algorithm();
}

void Sequencer::run(int numberOfCycles)
{
    stepsPerCycle = simParams->stepsPerCycle;
    if ( numberOfCycles ) this->numberOfCycles = numberOfCycles;
    else this->numberOfCycles = simParams->N; // / stepsPerCycle;
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
    int seq = 0;
    patch->positionsReady();
    suspend();
    DebugM(4,"Submit seq=" << seq << " Patch=" << patch->getPatchID() << "\n");
    reduction->submit(seq,REDUCTION_KINETIC_ENERGY,patch->calcKineticEnergy());
    // collection->submitPositions(seq,patch->atomIDList,patch->p);
    ++seq;
    for ( step = 0; step < numberOfCycles; ++step )
    {
	// DebugM(4,"Cycle #" << cycle << "\n");
        // for ( step = 0; step < stepsPerCycle; ++step )
        // {
            patch->addForceToMomentum(0.5*timestep);
            patch->addVelocityToPosition(timestep);
	    DebugM(4, patch->getPatchID()
		<< ": (" << cycle << "," << step << ") "
		<< "Sending positionsReady().\n");
            patch->positionsReady(!(step%stepsPerCycle));
            // patch->positionsReady(0);
	    DebugM(4, patch->getPatchID()
		<< ": (" << cycle << "," << step << ") "
		<< "Suspending.\n");
            suspend();
	    DebugM(4, patch->getPatchID()
		<< ": (" << cycle << "," << step << ") "
		<< "Awakened!\n");
            patch->addForceToMomentum(0.5*timestep);
	    DebugM(1,"Submit seq=" << seq << " Patch=" << patch->getPatchID() << "\n");
	    reduction->submit(seq, REDUCTION_KINETIC_ENERGY,
		patch->calcKineticEnergy());
	    // collection->submitPositions(seq,patch->atomIDList,patch->p);
	    ++seq;
        // }
    }
    DebugM(4, patch->getPatchID() << ": Exiting.\n");
    terminate();
}

