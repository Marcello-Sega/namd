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

static char ident[] = "@(#)$Header: /home/cvs/namd/cvsroot/namd2/src/Sequencer.C,v 1.1008 1997/02/14 05:53:04 jim Exp $";

#include "Node.h"
#include "SimParameters.h"
#include "Sequencer.h"
#include "HomePatch.h"
#include "ReductionMgr.h"
#include "CollectionMgr.h"

#define MIN_DEBUG_LEVEL 4
//#define DEBUGM
#include "Debug.h"

#define MIGRATION 1

Sequencer::Sequencer(HomePatch *p) :
	patch(p),
	simParams(Node::Object()->simParameters),
	reduction(ReductionMgr::Object()),
	collection(CollectionMgr::Object())
{
    reduction->Register(REDUCTION_KINETIC_ENERGY);
    threadStatus = NOTSUSPENDED;
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
    threadStatus = SUSPENDED;
    if ( numberOfCycles ) this->numberOfCycles = numberOfCycles;
    else this->numberOfCycles = simParams->N; // / stepsPerCycle;
    thread = CthCreate((CthVoidFn)&(threadRun),(void*)(this),0);
    CthSetStrategyDefault(thread);
    awaken();
}

void Sequencer::algorithm(void)
{
    const int numberOfCycles = this->numberOfCycles;
    const int stepsPerCycle = this->stepsPerCycle;
    const BigReal timestep = simParams->dt;
    int step, cycle=-1;	// cycle is unused!
    int seq = 0;
    // threadStatus = NOTSUSPENDED;
    patch->positionsReady();
    // if (threadStatus != AWAKENED)
//	{
	suspend();
//	}
    DebugM(4,"Submit seq=" << seq << " Patch=" << patch->getPatchID() << "\n");
    reduction->submit(seq,REDUCTION_KINETIC_ENERGY,patch->calcKineticEnergy());
    collection->submitPositions(seq,patch->atomIDList,patch->p);
    collection->submitVelocities(seq,patch->atomIDList,patch->v);
    ++seq;
    for ( step = 0; step < numberOfCycles; ++step )
    {
	// DebugM(4,"Cycle #" << step << "\n");
        // for ( step = 0; step < stepsPerCycle; ++step )
        // {
            patch->addForceToMomentum(0.5*timestep);
            patch->addVelocityToPosition(timestep);
	    DebugM(4, patch->getPatchID()
		<< ": (" << cycle << "," << step << ") "
		<< "Sending positionsReady().\n");
	    threadStatus = NOTSUSPENDED;
#if MIGRATION	// defined at the top of the file
            patch->positionsReady(!(step%stepsPerCycle));
#else
            patch->positionsReady(0);
#endif
            //Node::Object()->throwSequencer();
	    DebugM(4, patch->getPatchID()
		<< ": (" << cycle << "," << step << ") "
		<< "Suspending " << CthSelf() << " @" << CmiTimer() << "\n");
//	    if (threadStatus != AWAKENED)
//		{
		suspend();
//		}
	    DebugM(4, patch->getPatchID()
		<< ": (" << cycle << "," << step << ") "
		<< "Awakened!\n");
            patch->addForceToMomentum(0.5*timestep);
	    DebugM(4,"Submit seq=" << seq << " Patch=" << patch->getPatchID() << "\n");
	    reduction->submit(seq, REDUCTION_KINETIC_ENERGY,
		patch->calcKineticEnergy());
	    collection->submitPositions(seq,patch->atomIDList,patch->p);
	    collection->submitVelocities(seq,patch->atomIDList,patch->v);
	    ++seq;
        // }
    }
    DebugM(4, patch->getPatchID() << ": Exiting.\n");
    terminate();
}

void
Sequencer::terminate() {
  Node::messageHomeDone();
  CthFree(thread);
  CthSuspend();
}

/***************************************************************************
 * RCS INFORMATION:
 *
 *      $RCSfile: Sequencer.C,v $
 *      $Author: jim $  $Locker:  $             $State: Exp $
 *      $Revision: 1.1008 $     $Date: 1997/02/14 05:53:04 $
 *
 ***************************************************************************
 * REVISION HISTORY:
 *
 * $Log: Sequencer.C,v $
 * Revision 1.1008  1997/02/14 05:53:04  jim
 * Position and velocity output should now work.
 *
 * Revision 1.1007  1997/02/13 16:17:20  ari
 * Intermediate debuging commit - working to fix deep bug in migration?
 *
 ***************************************************************************/
