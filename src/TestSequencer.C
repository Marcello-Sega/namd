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

static char ident[] = "@(#)$Header: /home/cvs/namd/cvsroot/namd2/src/Attic/TestSequencer.C,v 1.2 1998/04/06 16:34:12 jim Exp $";

#include "Node.h"
#include "SimParameters.h"
#include "TestSequencer.h"
#include "HomePatch.h"
#include "ReductionMgr.h"
#include "CollectionMgr.h"
#include "BroadcastObject.h"
#include "Output.h"
#include "Controller.h"
#include "Broadcasts.h"
#include "Molecule.h"
#include "NamdOneTools.h"
#include "LdbCoordinator.h"
#include "Thread.h"

#define MIN_DEBUG_LEVEL 4
//#define DEBUGM
#include "Debug.h"

TestSequencer::TestSequencer(HomePatch *p) :
	Sequencer(p)
{
  ;
}

TestSequencer::~TestSequencer(void)
{
  ;
}

// Defines sequence of operations on a patch.  e.g. when
// to push out information for Compute objects to consume
// when to migrate atoms, when to add forces to velocity update.
void TestSequencer::algorithm(void)
{
    int &step = patch->flags.seq;
    step = simParams->firstTimestep;

    int &maxForceUsed = patch->flags.maxForceUsed;
    int &maxForceMerged = patch->flags.maxForceMerged;
    maxForceUsed = Results::normal;
    maxForceMerged = Results::normal;

    const int numberOfSteps = simParams->N;
    const int stepsPerCycle = simParams->stepsPerCycle;
    const BigReal timestep = simParams->dt;

    const int nonbondedFrequency = simParams->nonbondedFrequency;
    const BigReal nbondstep = timestep * nonbondedFrequency;
    int &doNonbonded = patch->flags.doNonbonded;
    doNonbonded = !(step%nonbondedFrequency);
    if ( nonbondedFrequency == 1 ) maxForceMerged = Results::nbond;
    if ( doNonbonded ) maxForceUsed = Results::nbond;

    // Do we do full electrostatics?
    const int dofull = ( simParams->fullDirectOn ||
			simParams->FMAOn || simParams->PMEOn );
    const int fullElectFrequency = simParams->fmaFrequency;
    const BigReal slowstep = timestep * fullElectFrequency;
    int &doFullElectrostatics = patch->flags.doFullElectrostatics;
    doFullElectrostatics = (dofull && !(step%fullElectFrequency));
    if ( dofull && (fullElectFrequency == 1) ) maxForceMerged = Results::slow;
    if ( doFullElectrostatics ) maxForceUsed = Results::slow;

    rattle1(0.);  // enforce rigid bond constraints on initial positions
    minimizationQuenchVelocity();
    runComputeObjects();
    addForceToMomentum(0.); // zero velocities of fixed atoms
    rattle2(timestep,step);  // enfore rigid bonds on initial velocities
    submitReductions(step);
    // submitCollections(step);
    rescaleVelocities(step);
    tcoupleVelocities(timestep,step);
    berendsenPressure(step);

    for ( ++step; step <= numberOfSteps; ++step )
    {
	langevinVelocities(0.5*timestep);

	/*
	addForceToMomentum(0.5*timestep);
	if (doNonbonded)
		addForceToMomentum(0.5*nbondstep,Results::nbond);
	if (doFullElectrostatics)
		addForceToMomentum(0.5*slowstep,Results::slow);
	*/

	minimizationMaximumMove(timestep);
	// addVelocityToPosition(timestep);
	translatePosition(0.2,0.2,0.2);
	rattle1(timestep);
	minimizationQuenchVelocity();

	doNonbonded = !(step%nonbondedFrequency);
	doFullElectrostatics = (dofull && !(step%fullElectFrequency));

        maxForceUsed = Results::normal;
	if ( doNonbonded ) maxForceUsed = Results::nbond;
	if ( doFullElectrostatics ) maxForceUsed = Results::slow;

	// Migrate Atoms on stepsPerCycle
	runComputeObjects(!(step%stepsPerCycle));

	/*
	addForceToMomentum(0.5*timestep);
	if (doNonbonded)
		addForceToMomentum(0.5*nbondstep,Results::nbond);
	if (doFullElectrostatics)
		addForceToMomentum(0.5*slowstep,Results::slow);
	*/

	langevinVelocities(0.5*timestep);

	rattle2(timestep,step);

	submitReductions(step);
	// submitCollections(step);
	rescaleVelocities(step);
	tcoupleVelocities(timestep,step);
	berendsenPressure(step);
#ifdef CYCLE_BARRIER
	int x;
	if (!((step+1) % stepsPerCycle))
	  x = broadcast->cycleBarrier.get(step);
#endif

	//
	// Trigger load balance stats collection
	//
	rebalanceLoad(step);
    }

    terminate();
}

void TestSequencer::translatePosition(BigReal dx, BigReal dy, BigReal dz) {
  Vector delta(dx,dy,dz);
  Position *p = patch->p.begin();
  Position *p_e = patch->p.end();
  for ( ; p != p_e; ++p ) *p += delta;
}


/***************************************************************************
 * RCS INFORMATION:
 *
 *      $RCSfile: TestSequencer.C,v $
 *      $Author: jim $  $Locker:  $             $State: Exp $
 *      $Revision: 1.2 $     $Date: 1998/04/06 16:34:12 $
 *
 ***************************************************************************
 * REVISION HISTORY:
 *
 * $Log: TestSequencer.C,v $
 * Revision 1.2  1998/04/06 16:34:12  jim
 * Added DPME (single processor only), test mode, and momenta printing.
 *
 * Revision 1.1  1998/03/31 04:55:50  jim
 * Added test mode, fixed errors in virial with full electrostatics.
 *
 *
 ***************************************************************************/
