/**
***  Copyright (c) 1995, 1996, 1997, 1998, 1999, 2000 by
***  The Board of Trustees of the University of Illinois.
***  All rights reserved.
**/

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
#include "Random.h"
#include "PatchMap.inl"

#define MIN_DEBUG_LEVEL 4
//#define DEBUGM
#include "Debug.h"

Sequencer::Sequencer(HomePatch *p) :
	simParams(Node::Object()->simParameters),
	patch(p),
	collection(CollectionMgr::Object())
{
    broadcast = new ControllerBroadcasts;
    reduction = ReductionMgr::Object()->willSubmit(REDUCTIONS_BASIC);
    ldbCoordinator = (LdbCoordinator::Object());
    random = new Random(simParams->randomSeed);
    random->split(patch->getPatchID()+1,PatchMap::Object()->numPatches()+1);

    rescaleVelocities_numTemps = 0;
}

Sequencer::~Sequencer(void)
{
    delete broadcast;
    delete reduction;
    delete random;
}

// Invoked by thread
void Sequencer::threadRun(Sequencer* arg)
{
    arg->algorithm();
}

// Invoked by Node::run() via HomePatch::runSequencer()
void Sequencer::run(void)
{
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
  int scriptTask = 1;
  int scriptSeq = 0;
  while ( (! simParams->tclOn) ||
		(scriptTask = broadcast->scriptBarrier.get(scriptSeq++)) ) {
    switch ( scriptTask ) {
      case 1:
	break;
      case 2:
	collection->submitPositions(0,patch->atomIDList,patch->p,
					patch->lattice,patch->t,2);
	continue;
      case 3:
	collection->submitVelocities(0,patch->atomIDList,patch->v);
	continue;
    }

    int &step = patch->flags.step;
    step = simParams->firstTimestep;

    int &maxForceUsed = patch->flags.maxForceUsed;
    int &maxForceMerged = patch->flags.maxForceMerged;
    maxForceUsed = Results::normal;
    maxForceMerged = Results::normal;

    const int numberOfSteps = simParams->N;
    const int stepsPerCycle = simParams->stepsPerCycle;
    const BigReal timestep = simParams->dt;

    const int nonbondedFrequency = simParams->nonbondedFrequency;
    slowFreq = nonbondedFrequency;
    const BigReal nbondstep = timestep * nonbondedFrequency;
    int &doNonbonded = patch->flags.doNonbonded;
    doNonbonded = !(step%nonbondedFrequency);
    if ( nonbondedFrequency == 1 ) maxForceMerged = Results::nbond;
    if ( doNonbonded ) maxForceUsed = Results::nbond;

    // Do we do full electrostatics?
    const int dofull = ( simParams->fullDirectOn ||
			simParams->FMAOn || simParams->PMEOn );
    const int fullElectFrequency = simParams->fullElectFrequency;
    if ( dofull ) slowFreq = fullElectFrequency;
    const BigReal slowstep = timestep * fullElectFrequency;
    int &doFullElectrostatics = patch->flags.doFullElectrostatics;
    doFullElectrostatics = (dofull && !(step%fullElectFrequency));
    if ( dofull && (fullElectFrequency == 1) && !(simParams->mollyOn) )
					maxForceMerged = Results::slow;
    if ( doFullElectrostatics ) maxForceUsed = Results::slow;
    int &doMolly = patch->flags.doMolly;
    doMolly = simParams->mollyOn && doFullElectrostatics;

    rattle1(0.);  // enforce rigid bond constraints on initial positions
    minimizationQuenchVelocity();
    runComputeObjects();
    addForceToMomentum(0.); // zero velocities of fixed atoms
    reassignVelocities(step);
    rattle2(timestep,step);  // enforce rigid bonds on initial velocities
    submitReductions(step);
    //submitCollections(step);
    //rescaleVelocities(step);
    //tcoupleVelocities(timestep,step);
    //berendsenPressure(step);

    for ( ++step; step <= numberOfSteps; ++step )
    {
	rescaleVelocities(step);
	tcoupleVelocities(timestep,step);
	berendsenPressure(step);
	// langevinVelocities(0.5*timestep);

	addForceToMomentum(0.5*timestep);
	if (doNonbonded)
		addForceToMomentum(0.5*nbondstep,Results::nbond);
	if (doFullElectrostatics)
		addForceToMomentum(0.5*slowstep,Results::slow);
	langevinVelocitiesBBK1(timestep);

	maximumMove(timestep);
	addVelocityToPosition(0.5*timestep);
	langevinPiston(step);
	addVelocityToPosition(0.5*timestep);
	rattle1(timestep);
	minimizationQuenchVelocity();

	doNonbonded = !(step%nonbondedFrequency);
	doFullElectrostatics = (dofull && !(step%fullElectFrequency));
	doMolly = simParams->mollyOn && doFullElectrostatics;

        maxForceUsed = Results::normal;
	if ( doNonbonded ) maxForceUsed = Results::nbond;
	if ( doFullElectrostatics ) maxForceUsed = Results::slow;

	// Migrate Atoms on stepsPerCycle
	runComputeObjects(!(step%stepsPerCycle));

	langevinVelocitiesBBK2(timestep);
	addForceToMomentum(0.5*timestep);
	if (doNonbonded)
		addForceToMomentum(0.5*nbondstep,Results::nbond);
	if (doFullElectrostatics)
		addForceToMomentum(0.5*slowstep,Results::slow);

	// langevinVelocities(0.5*timestep);
	reassignVelocities(step);

	rattle2(timestep,step);

	submitReductions(step);
	submitCollections(step);
	//rescaleVelocities(step);
	//tcoupleVelocities(timestep,step);
	//berendsenPressure(step);
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

    if (! simParams->tclOn) break;
  }

  submitCollections(0);
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
      patch->v[i] += f2 * random->gaussian_vector();
    }
  }
}

void Sequencer::langevinVelocitiesBBK1(BigReal dt_fs)
{
  if ( simParams->langevinOn )
  {
    Molecule *molecule = Node::Object()->molecule;
    BigReal dt = dt_fs * 0.001;  // convert to ps
    BigReal kbT = BOLTZMAN*(simParams->langevinTemp);
    for ( int i = 0; i < patch->numAtoms; ++i )
    {
      int aid = patch->atomIDList[i];
      BigReal dt_gamma = dt * molecule->langevin_param(aid);

      patch->v[i] += random->gaussian_vector() *
             sqrt( dt_gamma * kbT / patch->a[i].mass );
      patch->v[i] /= ( 1. + 0.5 * dt_gamma );
    }
  }
}


void Sequencer::langevinVelocitiesBBK2(BigReal dt_fs)
{
  if ( simParams->langevinOn )
  {
    Molecule *molecule = Node::Object()->molecule;
    BigReal dt = dt_fs * 0.001;  // convert to ps
    BigReal kbT = BOLTZMAN*(simParams->langevinTemp);
    for ( int i = 0; i < patch->numAtoms; ++i )
    {
      int aid = patch->atomIDList[i];
      BigReal dt_gamma = dt * molecule->langevin_param(aid);

      patch->v[i] *= ( 1. - 0.5 * dt_gamma );
      patch->v[i] += random->gaussian_vector() *
             sqrt( dt_gamma * kbT / patch->a[i].mass );
    }
  }
}

void Sequencer::berendsenPressure(int step)
{
  const int freq = simParams->berendsenPressureFreq;
  if ( simParams->berendsenPressureOn && !((step-1)%freq) )
  {
   Tensor factor = broadcast->positionRescaleFactor.get(step);
   patch->lattice.rescale(factor);
   if ( simParams->useGroupPressure )
   {
    int hgs;
    for ( int i = 0; i < patch->numAtoms; i += hgs ) {
      hgs = patch->a[i].hydrogenGroupSize;
      int j;
      BigReal m_cm = 0;
      Position x_cm;
      for ( j = i; j < (i+hgs); ++j ) {
        m_cm += patch->a[j].mass;
        x_cm += patch->a[j].mass * patch->p[j];
      }
      x_cm /= m_cm;
      Position new_x_cm = x_cm;
      patch->lattice.rescale(new_x_cm,factor);
      Position delta_x_cm = new_x_cm - x_cm;
      for ( j = i; j < (i+hgs); ++j ) {
        patch->p[j] += delta_x_cm;
      }
    }
   }
   else
   {
    for ( int i = 0; i < patch->numAtoms; ++i )
    {
      patch->lattice.rescale(patch->p[i],factor);
    }
   }
  }
}

void Sequencer::langevinPiston(int step)
{
  if ( simParams->langevinPistonOn && ! ( (step-1-slowFreq/2) % slowFreq ) )
  {
   Tensor factor = broadcast->positionRescaleFactor.get(step);
   // JCP FIX THIS!!!
   Vector velFactor(1/factor.xx,1/factor.yy,1/factor.zz);
   patch->lattice.rescale(factor);
   if ( simParams->useGroupPressure )
   {
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
      Position new_x_cm = x_cm;
      patch->lattice.rescale(new_x_cm,factor);
      Position delta_x_cm = new_x_cm - x_cm;
      v_cm /= m_cm;
      Velocity delta_v_cm;
      delta_v_cm.x = ( velFactor.x - 1 ) * v_cm.x;
      delta_v_cm.y = ( velFactor.y - 1 ) * v_cm.y;
      delta_v_cm.z = ( velFactor.z - 1 ) * v_cm.z;
      for ( j = i; j < (i+hgs); ++j ) {
        patch->p[j] += delta_x_cm;
        patch->v[j] += delta_v_cm;
      }
    }
   }
   else
   {
    for ( int i = 0; i < patch->numAtoms; ++i )
    {
      patch->lattice.rescale(patch->p[i],factor);
      patch->v[i].x *= velFactor.x;
      patch->v[i].y *= velFactor.y;
      patch->v[i].z *= velFactor.z;
    }
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
      rescaleVelocities_numTemps = 0;
    }
  }
}

void Sequencer::reassignVelocities(int step)
{
  const int reassignFreq = simParams->reassignFreq;
  if ( ( reassignFreq > 0 ) && ! ( step % reassignFreq ) ) {
    BigReal newTemp = simParams->reassignTemp;
    newTemp += ( step / reassignFreq ) * simParams->reassignIncr;
    if ( simParams->reassignIncr > 0.0 ) {
      if ( newTemp > simParams->reassignHold && simParams->reassignHold > 0.0 )
        newTemp = simParams->reassignHold;
    } else {
      if ( newTemp < simParams->reassignHold )
        newTemp = simParams->reassignHold;
    }
    BigReal kbT = BOLTZMAN * newTemp;
    for ( int i = 0; i < patch->numAtoms; ++i )
    {
      patch->v[i] = ( ( patch->a[i].flags & ATOM_FIXED ) ? Vector(0,0,0) :
        sqrt( kbT / patch->a[i].mass ) * random->gaussian_vector() );
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
    Vector virial(0.,0.,0.);
    patch->rattle2(dt, &virial);
    reduction->item(REDUCTION_VIRIAL_NORMAL_X) += virial.x;
    reduction->item(REDUCTION_VIRIAL_NORMAL_Y) += virial.y;
    reduction->item(REDUCTION_VIRIAL_NORMAL_Z) += virial.z;
    reduction->item(REDUCTION_ALT_VIRIAL_NORMAL_X) += virial.x;
    reduction->item(REDUCTION_ALT_VIRIAL_NORMAL_Y) += virial.y;
    reduction->item(REDUCTION_ALT_VIRIAL_NORMAL_Z) += virial.z;
    reduction->item(REDUCTION_INT_VIRIAL_NORMAL_X) += virial.x;
    reduction->item(REDUCTION_INT_VIRIAL_NORMAL_Y) += virial.y;
    reduction->item(REDUCTION_INT_VIRIAL_NORMAL_Z) += virial.z;
  }
}

void Sequencer::maximumMove(BigReal timestep)
{
  if ( simParams->maximumMove ) {
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
  } else {
    const BigReal dt = timestep / TIMEFACTOR;
    const BigReal maxvel = 10.0 / dt;
    const BigReal maxvel2 = maxvel * maxvel;
    int killme = 0;
    VelocityList::iterator v_i, v_e;
    v_i = patch->v.begin();  v_e = patch->v.end();
    for ( ; v_i != v_e; ++v_i ) {
      killme = killme || ( v_i->length2() > maxvel2 );
    }
    if ( killme ) {
      NAMD_die("Atoms moving too fast; simulation has become unstable.");
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
  reduction->item(REDUCTION_ATOM_CHECKSUM) += patch->getNumAtoms();
  reduction->item(REDUCTION_KINETIC_ENERGY) += patch->calcKineticEnergy();

  {
    Vector virial(0.,0.,0.);
    for ( int i = 0; i < patch->numAtoms; ++i ) {
      virial.x += ( patch->a[i].mass * patch->v[i].x * patch->v[i].x );
      virial.y += ( patch->a[i].mass * patch->v[i].y * patch->v[i].y );
      virial.z += ( patch->a[i].mass * patch->v[i].z * patch->v[i].z );
    }
    reduction->item(REDUCTION_VIRIAL_NORMAL_X) += virial.x;
    reduction->item(REDUCTION_VIRIAL_NORMAL_Y) += virial.y;
    reduction->item(REDUCTION_VIRIAL_NORMAL_Z) += virial.z;
  }
  {
    Vector altVirial(0.,0.,0.);
    for ( int i = 0; i < patch->numAtoms; ++i ) {
      altVirial.x += ( patch->f[Results::normal][i].x * patch->p[i].x );
      altVirial.y += ( patch->f[Results::normal][i].y * patch->p[i].y );
      altVirial.z += ( patch->f[Results::normal][i].z * patch->p[i].z );
      altVirial.x += ( patch->a[i].mass * patch->v[i].x * patch->v[i].x );
      altVirial.y += ( patch->a[i].mass * patch->v[i].y * patch->v[i].y );
      altVirial.z += ( patch->a[i].mass * patch->v[i].z * patch->v[i].z );
    }
    reduction->item(REDUCTION_ALT_VIRIAL_NORMAL_X) += altVirial.x;
    reduction->item(REDUCTION_ALT_VIRIAL_NORMAL_Y) += altVirial.y;
    reduction->item(REDUCTION_ALT_VIRIAL_NORMAL_Z) += altVirial.z;
  }
  {
    Vector altVirial(0.,0.,0.);
    for ( int i = 0; i < patch->numAtoms; ++i ) {
      altVirial.x += ( patch->f[Results::nbond][i].x * patch->p[i].x );
      altVirial.y += ( patch->f[Results::nbond][i].y * patch->p[i].y );
      altVirial.z += ( patch->f[Results::nbond][i].z * patch->p[i].z );
    }
    reduction->item(REDUCTION_ALT_VIRIAL_NBOND_X) += altVirial.x;
    reduction->item(REDUCTION_ALT_VIRIAL_NBOND_Y) += altVirial.y;
    reduction->item(REDUCTION_ALT_VIRIAL_NBOND_Z) += altVirial.z;
  }
  {
    Vector altVirial(0.,0.,0.);
    for ( int i = 0; i < patch->numAtoms; ++i ) {
      altVirial.x += ( patch->f[Results::slow][i].x * patch->p[i].x );
      altVirial.y += ( patch->f[Results::slow][i].y * patch->p[i].y );
      altVirial.z += ( patch->f[Results::slow][i].z * patch->p[i].z );
    }
    reduction->item(REDUCTION_ALT_VIRIAL_SLOW_X) += altVirial.x;
    reduction->item(REDUCTION_ALT_VIRIAL_SLOW_Y) += altVirial.y;
    reduction->item(REDUCTION_ALT_VIRIAL_SLOW_Z) += altVirial.z;
  }

  {
    BigReal intKineticEnergy = 0;
    Vector intVirialNormal(0.,0.,0.);
    Vector intVirialNbond(0.,0.,0.);
    Vector intVirialSlow(0.,0.,0.);

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
        intVirialNormal.x += patch->a[j].mass * (patch->v[j].x * dv.x);
        intVirialNormal.y += patch->a[j].mass * (patch->v[j].y * dv.y);
        intVirialNormal.z += patch->a[j].mass * (patch->v[j].z * dv.z);
        Vector dx = patch->p[j] - x_cm;
        intVirialNormal.x += patch->f[Results::normal][j].x * dx.x;
        intVirialNormal.y += patch->f[Results::normal][j].y * dx.y;
        intVirialNormal.z += patch->f[Results::normal][j].z * dx.z;
        intVirialNbond.x += patch->f[Results::nbond][j].x * dx.x;
        intVirialNbond.y += patch->f[Results::nbond][j].y * dx.y;
        intVirialNbond.z += patch->f[Results::nbond][j].z * dx.z;
        intVirialSlow.x += patch->f[Results::slow][j].x * dx.x;
        intVirialSlow.y += patch->f[Results::slow][j].y * dx.y;
        intVirialSlow.z += patch->f[Results::slow][j].z * dx.z;
      }
    }

    intKineticEnergy *= 0.5;

    reduction->item(REDUCTION_INT_KINETIC_ENERGY) += intKineticEnergy;
    reduction->item(REDUCTION_INT_VIRIAL_NORMAL_X) += intVirialNormal.x;
    reduction->item(REDUCTION_INT_VIRIAL_NORMAL_Y) += intVirialNormal.y;
    reduction->item(REDUCTION_INT_VIRIAL_NORMAL_Z) += intVirialNormal.z;
    reduction->item(REDUCTION_INT_VIRIAL_NBOND_X) += intVirialNbond.x;
    reduction->item(REDUCTION_INT_VIRIAL_NBOND_Y) += intVirialNbond.y;
    reduction->item(REDUCTION_INT_VIRIAL_NBOND_Z) += intVirialNbond.z;
    reduction->item(REDUCTION_INT_VIRIAL_SLOW_X) += intVirialSlow.x;
    reduction->item(REDUCTION_INT_VIRIAL_SLOW_Y) += intVirialSlow.y;
    reduction->item(REDUCTION_INT_VIRIAL_SLOW_Z) += intVirialSlow.z;
  }

  Vector momentum = patch->calcMomentum();
  reduction->item(REDUCTION_MOMENTUM_X) += momentum.x;
  reduction->item(REDUCTION_MOMENTUM_Y) += momentum.y;
  reduction->item(REDUCTION_MOMENTUM_Z) += momentum.z;

  Vector angularMomentum = patch->calcAngularMomentum();
  reduction->item(REDUCTION_ANGULAR_MOMENTUM_X) += angularMomentum.x;  
  reduction->item(REDUCTION_ANGULAR_MOMENTUM_Y) += angularMomentum.y;  
  reduction->item(REDUCTION_ANGULAR_MOMENTUM_Z) += angularMomentum.z;  

  reduction->submit();
}

void Sequencer::submitCollections(int step)
{
  if ( int prec = Output::coordinateNeeded(step) )
    collection->submitPositions(step,patch->atomIDList,patch->p,
				patch->lattice,patch->t,prec);
  if ( Output::velocityNeeded(step) )
    collection->submitVelocities(step,patch->atomIDList,patch->v);
}

void Sequencer::runComputeObjects(int migration)
{
  patch->positionsReady(migration);
  suspend(); // until all deposit boxes close
  if ( patch->flags.doMolly ) {
    Vector virial(0.,0.,0.);
    patch->mollyMollify(&virial);
    reduction->item(REDUCTION_VIRIAL_SLOW_X) += virial.x;
    reduction->item(REDUCTION_VIRIAL_SLOW_Y) += virial.y;
    reduction->item(REDUCTION_VIRIAL_SLOW_Z) += virial.z;
    reduction->item(REDUCTION_ALT_VIRIAL_SLOW_X) += virial.x;
    reduction->item(REDUCTION_ALT_VIRIAL_SLOW_Y) += virial.y;
    reduction->item(REDUCTION_ALT_VIRIAL_SLOW_Z) += virial.z;
    reduction->item(REDUCTION_INT_VIRIAL_SLOW_X) += virial.x;
    reduction->item(REDUCTION_INT_VIRIAL_SLOW_Y) += virial.y;
    reduction->item(REDUCTION_INT_VIRIAL_SLOW_Z) += virial.z;
  }
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
  CthFree(thread);
  CthSuspend();
}

