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

static char ident[] = "@(#)$Header: /home/cvs/namd/cvsroot/namd2/src/Sequencer.C,v 1.1047 1998/08/18 23:27:44 jim Exp $";

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
    reduction->Register(REDUCTION_INT_KINETIC_ENERGY);
    reduction->Register(REDUCTION_BC_ENERGY); // in case not used elsewhere
    reduction->Register(REDUCTION_SMD_ENERGY); // in case not used elsewhere
    if ( simParams->rigidBonds != RIGID_NONE ) {
	;// reduction->Register(REDUCTION_VIRIAL);
    }
    reduction->Register(REDUCTION_ALT_VIRIAL_NORMAL);
    reduction->Register(REDUCTION_ALT_VIRIAL_NBOND);
    reduction->Register(REDUCTION_ALT_VIRIAL_SLOW);
    reduction->Register(REDUCTION_INT_VIRIAL_NORMAL);
    reduction->Register(REDUCTION_INT_VIRIAL_NBOND);
    reduction->Register(REDUCTION_INT_VIRIAL_SLOW);
    reduction->Register(REDUCTION_MOMENTUM_X);
    reduction->Register(REDUCTION_MOMENTUM_Y);
    reduction->Register(REDUCTION_MOMENTUM_Z);
    reduction->Register(REDUCTION_ANGULAR_MOMENTUM_X);
    reduction->Register(REDUCTION_ANGULAR_MOMENTUM_Y);
    reduction->Register(REDUCTION_ANGULAR_MOMENTUM_Z);
    ldbCoordinator = (LdbCoordinator::Object());
}

Sequencer::~Sequencer(void)
{
    delete broadcast;

    reduction->unRegister(REDUCTION_KINETIC_ENERGY);
    reduction->unRegister(REDUCTION_INT_KINETIC_ENERGY);
    reduction->unRegister(REDUCTION_BC_ENERGY); // in case not used elsewhere
    reduction->unRegister(REDUCTION_SMD_ENERGY); // in case not used elsewhere
    if ( simParams->rigidBonds != RIGID_NONE ) {
	;// reduction->unRegister(REDUCTION_VIRIAL);
    }
    reduction->unRegister(REDUCTION_ALT_VIRIAL_NORMAL);
    reduction->unRegister(REDUCTION_ALT_VIRIAL_NBOND);
    reduction->unRegister(REDUCTION_ALT_VIRIAL_SLOW);
    reduction->unRegister(REDUCTION_INT_VIRIAL_NORMAL);
    reduction->unRegister(REDUCTION_INT_VIRIAL_NBOND);
    reduction->unRegister(REDUCTION_INT_VIRIAL_SLOW);
    reduction->unRegister(REDUCTION_MOMENTUM_X);
    reduction->unRegister(REDUCTION_MOMENTUM_Y);
    reduction->unRegister(REDUCTION_MOMENTUM_Z);
    reduction->unRegister(REDUCTION_ANGULAR_MOMENTUM_X);
    reduction->unRegister(REDUCTION_ANGULAR_MOMENTUM_Y);
    reduction->unRegister(REDUCTION_ANGULAR_MOMENTUM_Z);
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
    reassignVelocities(step);
    rattle2(timestep,step);  // enfore rigid bonds on initial velocities
    submitReductions(step);
    submitCollections(step);
    rescaleVelocities(step);
    tcoupleVelocities(timestep,step);
    berendsenPressure(step);

    for ( ++step; step <= numberOfSteps; ++step )
    {
	langevinVelocities(0.5*timestep);

	addForceToMomentum(0.5*timestep);
	if (doNonbonded)
		addForceToMomentum(0.5*nbondstep,Results::nbond);
	if (doFullElectrostatics)
		addForceToMomentum(0.5*slowstep,Results::slow);

	minimizationMaximumMove(timestep);
	addVelocityToPosition(0.5*timestep);
	langevinPiston(step);
	addVelocityToPosition(0.5*timestep);
	rattle1(timestep);
	minimizationQuenchVelocity();

	doNonbonded = !(step%nonbondedFrequency);
	doFullElectrostatics = (dofull && !(step%fullElectFrequency));

        maxForceUsed = Results::normal;
	if ( doNonbonded ) maxForceUsed = Results::nbond;
	if ( doFullElectrostatics ) maxForceUsed = Results::slow;

	// Migrate Atoms on stepsPerCycle
	runComputeObjects(!(step%stepsPerCycle));

	addForceToMomentum(0.5*timestep);
	if (doNonbonded)
		addForceToMomentum(0.5*nbondstep,Results::nbond);
	if (doFullElectrostatics)
		addForceToMomentum(0.5*slowstep,Results::slow);

	langevinVelocities(0.5*timestep);
	reassignVelocities(step);

	rattle2(timestep,step);

	submitReductions(step);
	submitCollections(step);
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

void Sequencer::langevinVelocities(BigReal dt_fs)
{
  if ( simParams->langevinOn )
  {
    Molecule *molecule = Node::Object()->molecule;
    BigReal dt = dt_fs * 0.001;  // convert to ps
    BigReal kbT = BOLTZMAN*(simParams->langevinTemp);
    for ( int i = 0; i < patch->numAtoms; ++i )
    {
      int aid = patch->atomIDList[i];
      BigReal f1 = exp( -1. * dt * molecule->langevin_param(aid) );
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

void Sequencer::langevinPiston(int step)
{
  if ( simParams->langevinPistonOn )
  {
    BigReal factor = broadcast->positionRescaleFactor.get(step);
    BigReal velFactor = 1 / factor;
    patch->lattice.rescale(factor);
    for ( int i = 0; i < patch->numAtoms; ++i )
    {
      patch->lattice.rescale(patch->p[i],factor);
      patch->v[i] *= velFactor;
    }
  }
}

void Sequencer::rescaleVelocities(int step)
{
  const int rescaleFreq = simParams->rescaleFreq;
  if ( rescaleFreq > 0 ) {
    ++rescaleVelocities_numTemps;
    if ( rescaleVelocities_numTemps == rescaleFreq ) {
      BigReal factor = broadcast->velocityRescaleFactor.get(step);
      for ( int i = 0; i < patch->numAtoms; ++i )
      {
        patch->v[i] *= factor;
      }
    }
  }
}

void Sequencer::reassignVelocities(int step)
{
  const int reassignFreq = simParams->reassignFreq;
  if ( ( reassignFreq > 0 ) && ! ( step % reassignFreq ) ) {
    BigReal newTemp = simParams->reassignTemp;
    newTemp += ( step / reassignFreq ) * simParams->reassignIncr;
    if ( newTemp < 0 ) newTemp = 0;
    BigReal kbT = BOLTZMAN * newTemp;
    for ( int i = 0; i < patch->numAtoms; ++i )
    {
      patch->v[i] = ( ( patch->a[i].flags & ATOM_FIXED ) ? Vector(0,0,0) :
        sqrt( kbT / patch->a[i].mass ) * gaussian_random_vector() );
    }
  }
}

void Sequencer::tcoupleVelocities(BigReal dt_fs, int step)
{
  if ( simParams->tCoupleOn )
  {
    BigReal coefficient = broadcast->tcoupleCoefficient.get(step);
    Molecule *molecule = Node::Object()->molecule;
    BigReal dt = dt_fs * 0.001;  // convert to ps
    coefficient *= dt;
    for ( int i = 0; i < patch->numAtoms; ++i )
    {
      int aid = patch->atomIDList[i];
      BigReal f1 = exp( coefficient * molecule->langevin_param(aid) );
      patch->v[i] *= f1;
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

void Sequencer::rattle1(BigReal dt)
{
  if ( simParams->rigidBonds != RIGID_NONE ) {
    patch->rattle1(dt);
  }
}

void Sequencer::rattle2(BigReal dt, int step)
{
  if ( simParams->rigidBonds != RIGID_NONE ) {
    BigReal virial = 0.;
    patch->rattle2(dt, &virial);
    // reduction->submit(step, REDUCTION_VIRIAL, virial);
  }
}

void Sequencer::minimizationMaximumMove(BigReal timestep)
{
  if ( simParams->minimizeOn ) {
    const BigReal dt = timestep / TIMEFACTOR;
    const BigReal maxvel = simParams->maximumMove / dt;
    const BigReal maxvel2 = maxvel * maxvel;
    VelocityList::iterator v_i, v_e;
    v_i = patch->v.begin();  v_e = patch->v.end();
    for ( ; v_i != v_e; ++v_i ) {
      if ( v_i->length2() > maxvel2 ) {
	*v_i *= ( maxvel / v_i->length() );
      }
    }
  }
}

void Sequencer::minimizationQuenchVelocity(void)
{
  if ( simParams->minimizeOn ) {
    VelocityList::iterator v_i, v_e;
    v_i = patch->v.begin();  v_e = patch->v.end();
    for ( ; v_i != v_e; ++v_i ) {
      *v_i = 0.;
    }
  }
}

void Sequencer::submitReductions(int step)
{
  reduction->submit(step,REDUCTION_KINETIC_ENERGY,patch->calcKineticEnergy());
  reduction->submit(step,REDUCTION_BC_ENERGY,0.);
  reduction->submit(step,REDUCTION_SMD_ENERGY,0.);  

  {
    BigReal altVirial = 0.;
    for ( int i = 0; i < patch->numAtoms; ++i ) {
      altVirial += ( patch->f[Results::normal][i] * patch->p[i] );
    }
    reduction->submit(step,REDUCTION_ALT_VIRIAL_NORMAL,altVirial);
  }
  {
    BigReal altVirial = 0.;
    for ( int i = 0; i < patch->numAtoms; ++i ) {
      altVirial += ( patch->f[Results::nbond][i] * patch->p[i] );
    }
    reduction->submit(step,REDUCTION_ALT_VIRIAL_NBOND,altVirial);
  }
  {
    BigReal altVirial = 0.;
    for ( int i = 0; i < patch->numAtoms; ++i ) {
      altVirial += ( patch->f[Results::slow][i] * patch->p[i] );
    }
    reduction->submit(step,REDUCTION_ALT_VIRIAL_SLOW,altVirial);
  }

  {
    BigReal intKineticEnergy = 0;
    BigReal intVirialNormal = 0;
    BigReal intVirialNbond = 0;
    BigReal intVirialSlow = 0;

    int hgs;
    for ( int i = 0; i < patch->numAtoms; i += hgs ) {
      hgs = patch->a[i].hydrogenGroupSize;
      int j;
      BigReal m_cm = 0;
      Position x_cm;
      Velocity v_cm;
      for ( j = i; j < (i+hgs); ++j ) {
        m_cm += patch->a[j].mass;
        x_cm += patch->a[j].mass * patch->p[j];
        v_cm += patch->a[j].mass * patch->v[j];
      }
      x_cm /= m_cm;
      v_cm /= m_cm;
      for ( j = i; j < (i+hgs); ++j ) {
        Vector dv = patch->v[j] - v_cm;
        intKineticEnergy += patch->a[j].mass * (patch->v[j] * dv);
        Vector dx = patch->p[j] - x_cm;
        intVirialNormal += patch->f[Results::normal][j] * dx;
        intVirialNbond += patch->f[Results::nbond][j] * dx;
        intVirialSlow += patch->f[Results::slow][j] * dx;
      }
    }

    intKineticEnergy *= 0.5;

    reduction->submit(step,REDUCTION_INT_KINETIC_ENERGY,intKineticEnergy);
    reduction->submit(step,REDUCTION_INT_VIRIAL_NORMAL,intVirialNormal);
    reduction->submit(step,REDUCTION_INT_VIRIAL_NBOND,intVirialNbond);
    reduction->submit(step,REDUCTION_INT_VIRIAL_SLOW,intVirialSlow);
  }

  Vector momentum = patch->calcMomentum();
  reduction->submit(step,REDUCTION_MOMENTUM_X,momentum.x);  
  reduction->submit(step,REDUCTION_MOMENTUM_Y,momentum.y);  
  reduction->submit(step,REDUCTION_MOMENTUM_Z,momentum.z);  

  Vector angularMomentum = patch->calcAngularMomentum();
  reduction->submit(step,REDUCTION_ANGULAR_MOMENTUM_X,angularMomentum.x);  
  reduction->submit(step,REDUCTION_ANGULAR_MOMENTUM_Y,angularMomentum.y);  
  reduction->submit(step,REDUCTION_ANGULAR_MOMENTUM_Z,angularMomentum.z);  

}

void Sequencer::submitCollections(int step)
{
  if ( Output::coordinateNeeded(step) )
    collection->submitPositions(step,patch->atomIDList,patch->p,
					patch->lattice,patch->t);
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
  if ( ldbCoordinator->balanceNow(timestep) )
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
 *      $Author: jim $  $Locker:  $             $State: Exp $
 *      $Revision: 1.1047 $     $Date: 1998/08/18 23:27:44 $
 *
 ***************************************************************************
 * REVISION HISTORY:
 *
 * $Log: Sequencer.C,v $
 * Revision 1.1047  1998/08/18 23:27:44  jim
 * First implementation of constant pressure.
 * Isotropic only, incompatible with multiple timestepping or SHAKE.
 *
 * Revision 1.1046  1998/08/11 16:30:31  jim
 * Modified output from periodic boundary simulations to return atoms to
 * internally consistent coordinates.  We store the transformations which
 * were performed and undo them at the end.  It might be better to do this
 * by always keeping the original coordinates and only doing the transform
 * for the nonbonded terms but this works for now.
 *
 * Revision 1.1045  1998/08/03 15:31:20  jim
 * Added temperature reassignment.
 *
 * Revision 1.1044  1998/08/02 21:26:40  jim
 * Altered velocity rescaling to use averaged temperature.
 *
 * Revision 1.1043  1998/06/18 14:48:05  jim
 * Split virial into NORMAL, NBOND, and SLOW parts to match force classes.
 *
 * Revision 1.1042  1998/04/06 16:34:09  jim
 * Added DPME (single processor only), test mode, and momenta printing.
 *
 * Revision 1.1041  1998/03/06 20:55:26  jim
 * Added temperature coupling.
 *
 * Revision 1.1040  1998/03/06 10:25:27  jim
 * Added very basic minimizer.
 *
 * Revision 1.1039  1998/02/18 19:13:59  jim
 * Fixed Langevin dynamics, undoing changes from yesterday.
 *
 * Revision 1.1038  1998/02/18 05:38:32  jim
 * RigidBonds mainly finished.  Now temperature is correct and a form
 * of Langevin dynamics works with constraints.
 *
 * Revision 1.1037  1998/02/17 06:39:23  jim
 * SHAKE/RATTLE (rigidBonds) appears to work!!!  Still needs langevin,
 * proper startup, and degree of freedom tracking.
 *
 * Revision 1.1036  1998/01/05 20:26:49  sergei
 * added reduction->(un)Register(REDUCTION_SMD_ENERGY) to (con/de)structor
 * added reduction->submit(step,REDUCTION_SMD_ENERGY,0.) to submitReductions()
 *
 * Revision 1.1035  1997/12/22 21:29:28  jim
 * Proxies no longer send empty arrays back to HomePatch.  Requires some new
 * flags to be set correctly in Sequencer in order to work.  These are:
 *   maxForceMerged - this and faster are added into Results::normal array
 *   maxForceUsed - all forces slower than this are discarded (assumed zero)
 * Generally maxForceMerged doesn't change but maxForceUsed depends on timestep.
 *
 * Revision 1.1034  1997/09/19 09:39:06  jim
 * Small tweaks for fixed atoms.
 *
 * Revision 1.1033  1997/09/19 08:55:36  jim
 * Added rudimentary but relatively efficient fixed atoms.  New options
 * are fixedatoms, fixedatomsfile, and fixedatomscol (nonzero means fixed).
 * Energies will be affected, although this can be fixed with a little work.
 *
 * Revision 1.1032  1997/08/22 19:27:37  brunner
 * Added cycle barrier, enabled by compiling with -DCYCLE_BARRIER
 *
 * Revision 1.1031  1997/04/16 22:12:20  brunner
 * Fixed an LdbCoordinator bug, and cleaned up timing and Ldb output some.
 *
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
