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

static char ident[] = "@(#)$Header: /home/cvs/namd/cvsroot/namd2/src/Sequencer.C,v 1.1019 1997/03/19 22:44:24 jim Exp $";

#include "Node.h"
#include "SimParameters.h"
#include "Sequencer.h"
#include "HomePatch.h"
#include "ReductionMgr.h"
#include "CollectionMgr.h"
#include "BroadcastObject.h"
#include "Output.h"
#include "Controller.h"

#define MIN_DEBUG_LEVEL 3
//#define DEBUGM
#include "Debug.h"

Sequencer::Sequencer(HomePatch *p) :
	patch(p),
	simParams(Node::Object()->simParameters),
	reduction(ReductionMgr::Object()),
	collection(CollectionMgr::Object())
{
    if ( simParams->rescaleFreq > 0 )
    {
	velocityRescaleFactor = new
		SimpleBroadcastObject<BigReal>(velocityRescaleFactorTag);
    }
    else velocityRescaleFactor = 0;

    reduction->Register(REDUCTION_KINETIC_ENERGY);
    reduction->Register(REDUCTION_BC_ENERGY); // in case not used elsewhere
}

Sequencer::~Sequencer(void)
{
    delete velocityRescaleFactor;

    reduction->unRegister(REDUCTION_KINETIC_ENERGY);
    reduction->unRegister(REDUCTION_BC_ENERGY); // in case not used elsewhere
}

// Invoked by thread
void Sequencer::threadRun(Sequencer* arg)
{
    arg->algorithm();
}

// Invoked by Node::run() via HomePatch::runSequencer()
void Sequencer::run(int numberOfCycles)
{
    this->numberOfCycles = numberOfCycles;
    if ( numberOfCycles ) 
      NAMD_die("Sorry, Sequencer::run() does not support an argument.\n");

    // create a Thread and invoke it
    DebugM(4, "::run() - this = " << this << "\n" );
    thread = CthCreate((CthVoidFn)&(threadRun),(void*)(this),0);
    CthSetStrategyDefault(thread);
    awaken();
}

// Defines sequence of operations on a patch.  e.g. when
// to push out information for Compute objects to consume
// when to migrate atoms, when to add forces to velocity update.
void Sequencer::algorithm(void)
{
    int &step = patch->flags.seq;
    step = simParams->firstTimestep;

    const int numberOfSteps = simParams->N;
    const int stepsPerCycle = simParams->stepsPerCycle;
    const BigReal timestep = simParams->dt;

    // Do we do full electrostatics?
    const int dofull = ( simParams->fullDirectOn || simParams->FMAOn );
    const BigReal slowstep = timestep * stepsPerCycle;
    int &doFullElectrostatics = patch->flags.doFullElectrostatics;
    doFullElectrostatics = dofull;

    // Push out inital positions
    patch->positionsReady();
    suspend(); // until all deposit boxes close

    reduction->submit(step,REDUCTION_KINETIC_ENERGY,patch->calcKineticEnergy());
    reduction->submit(step,REDUCTION_BC_ENERGY,0.);
    submitCollections(step);
    rescaleVelocities(step);

    for ( ++step; step <= numberOfSteps; ++step )
    {
	patch->addForceToMomentum(0.5*timestep);
	if (dofull && !(step%stepsPerCycle))
		patch->addForceToMomentum(0.5*slowstep,Results::slow);
	patch->addVelocityToPosition(timestep);

	doFullElectrostatics = (dofull && !(step%stepsPerCycle));

	// Migrate Atoms on stepsPerCycle
	patch->positionsReady(!(step%stepsPerCycle));
	suspend(); // until all deposit boxes close

	patch->addForceToMomentum(0.5*timestep);
	if (dofull && !(step%stepsPerCycle))
		patch->addForceToMomentum(0.5*slowstep,Results::slow);

	// Pass up information from this Patch
	reduction->submit(step,REDUCTION_KINETIC_ENERGY,patch->calcKineticEnergy());
	reduction->submit(step,REDUCTION_BC_ENERGY,0.);
	submitCollections(step);
        rescaleVelocities(step);
    }

    terminate();
}

void Sequencer::rescaleVelocities(int step)
{
  const int rescaleFreq = simParams->rescaleFreq;
  if ( rescaleFreq > 0 && !(step%rescaleFreq) )
  {
    BigReal factor = velocityRescaleFactor->get(step);
    for ( int i = 0; i < patch->numAtoms; ++i )
    {
      patch->v[i] *= factor;
    }
  }
}

void Sequencer::submitCollections(int timestep)
{
  if ( Output::coordinateNeeded(timestep) )
    collection->submitPositions(timestep,patch->atomIDList,patch->p);
  if ( Output::velocityNeeded(timestep) )
    collection->submitVelocities(timestep,patch->atomIDList,patch->v);
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
 *      $Revision: 1.1019 $     $Date: 1997/03/19 22:44:24 $
 *
 ***************************************************************************
 * REVISION HISTORY:
 *
 * $Log: Sequencer.C,v $
 * Revision 1.1019  1997/03/19 22:44:24  jim
 * Revamped Controller/Sequencer, added velocity rescaling.
 *
 * Revision 1.1018  1997/03/19 11:54:55  ari
 * Add Broadcast mechanism.
 * Fixed RCS Log entries on files that did not have Log entries.
 * Added some register variables to Molecule and ComputeNonbondedExcl.C
 *
 * Revision 1.1017  1997/03/18 21:35:35  jim
 * Eliminated fake_seq.  Reductions now use Patch::flags.seq.
 *
 * Revision 1.1016  1997/03/18 18:09:15  jim
 * Revamped collection system to ensure ordering and eliminate
 * unnecessary collections.  Also reduced make dependencies.
 *
 * Revision 1.1015  1997/03/15 22:15:29  jim
 * Added ComputeCylindricalBC.  Doesn't break anything but untested and
 * cylinder is along x axis (will fix soon).
 *
 * Revision 1.1014  1997/03/14 21:40:15  ari
 * Reorganized startup to make possible inital load
 * balancing by changing methods in WorkDistrib.
 * Also made startup more transparent and easier
 * to modify.
 *
 * Revision 1.1013  1997/03/13 06:37:12  jim
 * Multiple time-stepping implemented, still needs proper splitting functions.
 *
 * Revision 1.1012  1997/03/11 07:30:20  jim
 * Once more change to support firstTimestep parameter.
 *
 * Revision 1.1011  1997/03/11 04:04:38  jim
 * Now properly handles firstTimeStep parameter.
 *
 * Revision 1.1010  1997/03/04 22:37:17  ari
 * Clean up of code.  Debug statements removal, dead code removal.
 * Minor fixes, output fixes.
 * Commented some code from the top->down.  Mainly reworked Namd, Node, main.
 *
 * Revision 1.1009  1997/02/28 04:47:13  jim
 * Full electrostatics now works with fulldirect on one node.
 *
 * Revision 1.1008  1997/02/14 05:53:04  jim
 * Position and velocity output should now work.
 *
 * Revision 1.1007  1997/02/13 16:17:20  ari
 * Intermediate debuging commit - working to fix deep bug in migration?
 *
 ***************************************************************************/
