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

static char ident[] = "@(#)$Header: /home/cvs/namd/cvsroot/namd2/src/Sequencer.C,v 1.1030 1997/04/10 09:14:12 ari Exp $";

#include "Node.h"
#include "SimParameters.h"
#include "Sequencer.h"
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

Sequencer::Sequencer(HomePatch *p) :
	patch(p),
	simParams(Node::Object()->simParameters),
	reduction(ReductionMgr::Object()),
	collection(CollectionMgr::Object())
{
    broadcast = new ControllerBroadcasts;

    reduction->Register(REDUCTION_KINETIC_ENERGY);
    reduction->Register(REDUCTION_BC_ENERGY); // in case not used elsewhere
    reduction->Register(REDUCTION_ALT_VIRIAL);
    ldbCoordinator = (LdbCoordinator::Object());
}

Sequencer::~Sequencer(void)
{
    delete broadcast;

    reduction->unRegister(REDUCTION_KINETIC_ENERGY);
    reduction->unRegister(REDUCTION_BC_ENERGY); // in case not used elsewhere
    reduction->unRegister(REDUCTION_ALT_VIRIAL);
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
    thread = CthCreate((CthVoidFn)&(threadRun),(void*)(this),SEQ_STK_SZ);
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
    doFullElectrostatics = (dofull && !(step%stepsPerCycle));

    const int nonbondedFrequency = simParams->nonbondedFrequency;
    const BigReal nbondstep = timestep * nonbondedFrequency;
    int &doNonbonded = patch->flags.doNonbonded;
    doNonbonded = !(step%nonbondedFrequency);

    runComputeObjects();
    submitReductions(step);
    submitCollections(step);
    rescaleVelocities(step);
    berendsenPressure(step);
    langevinVelocities(step);

    for ( ++step; step <= numberOfSteps; ++step )
    {
	addForceToMomentum(0.5*timestep);
	if (doNonbonded)
		addForceToMomentum(0.5*nbondstep,Results::nbond);
	if (doFullElectrostatics)
		addForceToMomentum(0.5*slowstep,Results::slow);

	addVelocityToPosition(timestep);

	doNonbonded = !(step%nonbondedFrequency);
	doFullElectrostatics = (dofull && !(step%stepsPerCycle));

	// Migrate Atoms on stepsPerCycle
	runComputeObjects(!(step%stepsPerCycle));

	addForceToMomentum(0.5*timestep);
	if (doNonbonded)
		addForceToMomentum(0.5*nbondstep,Results::nbond);
	if (doFullElectrostatics)
		addForceToMomentum(0.5*slowstep,Results::slow);

	submitReductions(step);
	submitCollections(step);
	rescaleVelocities(step);
	berendsenPressure(step);
	langevinVelocities(step);
	//
	// Trigger load balance stats collection
	//
	rebalanceLoad(step);
    }

    terminate();
}

void Sequencer::langevinVelocities(int step)
{
  if ( simParams->langevinOn )
  {
    Molecule *molecule = Node::Object()->molecule;
    BigReal dt = simParams->dt * 0.001;  // convert to ps
    BigReal kbT = BOLTZMAN*(simParams->langevinTemp);
    for ( int i = 0; i < patch->numAtoms; ++i )
    {
      int aid = patch->atomIDList[i];
      BigReal f1 = exp( -1. * dt * molecule->langevin_param(aid) );
      DebugM(1,"At step " << step << " langevin decay is " << f1 << "\n");
      BigReal f2 = sqrt( ( 1. - f1*f1 ) * kbT / patch->a[i].mass );

      patch->v[i] *= f1;
      patch->v[i] += f2 * gaussian_random_vector();
    }
  }
}

void Sequencer::berendsenPressure(int step)
{
  const int freq = simParams->berendsenPressureFreq;
  if ( simParams->berendsenPressureOn && !(step%freq) )
  {
    BigReal factor = broadcast->positionRescaleFactor.get(step);
    patch->lattice.rescale(factor);
    for ( int i = 0; i < patch->numAtoms; ++i )
    {
      patch->lattice.rescale(patch->p[i],factor);
    }
  }
}

void Sequencer::rescaleVelocities(int step)
{
  const int rescaleFreq = simParams->rescaleFreq;
  if ( rescaleFreq > 0 && !(step%rescaleFreq) )
  {
    BigReal factor = broadcast->velocityRescaleFactor.get(step);
    for ( int i = 0; i < patch->numAtoms; ++i )
    {
      patch->v[i] *= factor;
    }
  }
}

void Sequencer::addForceToMomentum(BigReal dt, const int ftag)
{
  patch->addForceToMomentum(dt,ftag);
}

void Sequencer::addVelocityToPosition(BigReal dt)
{
  patch->addVelocityToPosition(dt);
}

void Sequencer::submitReductions(int step)
{
  reduction->submit(step,REDUCTION_KINETIC_ENERGY,patch->calcKineticEnergy());
  reduction->submit(step,REDUCTION_BC_ENERGY,0.);

  BigReal altVirial = 0.;
  for ( int i = 0; i < patch->numAtoms; ++i )
  {
    for ( int j = 0; j < Results::maxNumForces; ++j )
    {
      altVirial += ( patch->f[j][i] * patch->p[i] );
    }
  }
  reduction->submit(step,REDUCTION_ALT_VIRIAL,altVirial);
}

void Sequencer::submitCollections(int step)
{
  if ( Output::coordinateNeeded(step) )
    collection->submitPositions(step,patch->atomIDList,patch->p);
  if ( Output::velocityNeeded(step) )
    collection->submitVelocities(step,patch->atomIDList,patch->v);
}

void Sequencer::runComputeObjects(int migration)
{
  patch->positionsReady(migration);
  suspend(); // until all deposit boxes close
}

void Sequencer::rebalanceLoad(int timestep)
{
  if ( (simParams->ldbStrategy != LDBSTRAT_NONE)
       && ((timestep % simParams->ldbStepsPerCycle) == 0))
  {
    DebugM(4, "Running ldb!\n");
    patch->submitLoadStats(timestep);
    ldbCoordinator->rebalance(this,patch->getPatchID());
  }
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
 *      $Author: ari $  $Locker:  $             $State: Exp $
 *      $Revision: 1.1030 $     $Date: 1997/04/10 09:14:12 $
 *
 ***************************************************************************
 * REVISION HISTORY:
 *
 * $Log: Sequencer.C,v $
 * Revision 1.1030  1997/04/10 09:14:12  ari
 * Final debugging for compute migration / proxy creation for load balancing.
 * Lots of debug code added, mostly turned off now.
 * Fixed bug in PositionBox when Patch had no dependencies.
 * Eliminated use of cout and misuse of iout in numerous places.
 *                                            Ari & Jim
 *
 * Revision 1.1029  1997/04/08 21:08:47  jim
 * Contant pressure now correct on multiple nodes, should work with MTS.
 *
 * Revision 1.1028  1997/04/08 07:09:01  ari
 * Modification for dynamic loadbalancing - moving computes
 * Still bug in new computes or usage of proxies/homepatches.
 * Works if ldbStrategy is none as before.
 *
 * Revision 1.1027  1997/04/01 23:20:17  brunner
 * Collection on node 0 added
 *
 * Revision 1.1026  1997/03/27 20:25:51  brunner
 * Changes for LdbCoordinator, the load balance control BOC
 *
 * Revision 1.1025  1997/03/27 08:04:23  jim
 * Reworked Lattice to keep center of cell fixed during rescaling.
 *
 * Revision 1.1024  1997/03/27 03:16:56  jim
 * Added code to check virial calculation, fixed problems with DPMTA and PBC's.
 *
 * Revision 1.1023  1997/03/25 23:01:02  jim
 * Added nonbondedFrequency parameter and multiple time-stepping
 *
 * Revision 1.1022  1997/03/25 04:04:57  jim
 * Simplified algorithm a bit.
 *
 * Revision 1.1021  1997/03/24 01:44:01  jim
 * Added Langevin dynamics.
 *
 * Revision 1.1020  1997/03/21 23:05:41  jim
 * Added Berendsen's pressure coupling method, won't work with MTS yet.
 *
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
