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
	collection(CollectionMgr::Object()),
	ldbSteps(0)
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
  int scriptTask = SCRIPT_RUN;
  int scriptSeq = 0;
  while ( (! simParams->tclOn) ||
    (scriptTask = broadcast->scriptBarrier.get(scriptSeq++)) != SCRIPT_END ) {
    switch ( scriptTask ) {
      case SCRIPT_RUN:
	break;
      case SCRIPT_OUTPUT:
	submitCollections(FILE_OUTPUT);
	continue;
      case SCRIPT_MEASURE:
	submitCollections(EVAL_MEASURE);
	continue;
      case SCRIPT_REINITVELS:
	reinitVelocities();
	continue;
      case SCRIPT_CHECKPOINT:
        patch->checkpoint();
	continue;
      case SCRIPT_REVERT:
        patch->revert();
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
    rebalanceLoad(step);

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
        cycleBarrier(!((step+1) % stepsPerCycle),step);
	rebalanceLoad(step);
    }

    if (! simParams->tclOn) break;
  }

  // only reach here on SCRIPT_END or no script
  submitCollections(END_OF_RUN);
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

void Sequencer::reinitVelocities(void)
{
  BigReal newTemp = simParams->initialTemp;
  BigReal kbT = BOLTZMAN * newTemp;
  for ( int i = 0; i < patch->numAtoms; ++i )
  {
    patch->v[i] = ( ( patch->a[i].flags & ATOM_FIXED ) ? Vector(0,0,0) :
      sqrt( kbT / patch->a[i].mass ) * random->gaussian_vector() );
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
    Tensor virial;
    patch->rattle2(dt, &virial);
    ADD_TENSOR_OBJECT(reduction,REDUCTION_VIRIAL_NORMAL,virial);
    // we need to add to alt and int virial because not included in forces
    ADD_TENSOR_OBJECT(reduction,REDUCTION_ALT_VIRIAL_NORMAL,virial);
    ADD_TENSOR_OBJECT(reduction,REDUCTION_INT_VIRIAL_NORMAL,virial);
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
      iout << iERROR << 
        "Atoms moving too fast; simulation has become unstable.\n" << endi;
      Node::Object()->enableHaltBarrier();
      terminate();
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
    Tensor virial;
    for ( int i = 0; i < patch->numAtoms; ++i ) {
      virial += ( patch->a[i].mass * outer(patch->v[i],patch->v[i]) );
    }
    ADD_TENSOR_OBJECT(reduction,REDUCTION_VIRIAL_NORMAL,virial);
    ADD_TENSOR_OBJECT(reduction,REDUCTION_ALT_VIRIAL_NORMAL,virial);
  }
  {
    Tensor altVirial;
    for ( int i = 0; i < patch->numAtoms; ++i ) {
      altVirial += outer(patch->f[Results::normal][i],patch->p[i]);
    }
    ADD_TENSOR_OBJECT(reduction,REDUCTION_ALT_VIRIAL_NORMAL,altVirial);
  }
  {
    Tensor altVirial;
    for ( int i = 0; i < patch->numAtoms; ++i ) {
      altVirial += outer(patch->f[Results::nbond][i],patch->p[i]);
    }
    ADD_TENSOR_OBJECT(reduction,REDUCTION_ALT_VIRIAL_NBOND,altVirial);
  }
  {
    Tensor altVirial;
    for ( int i = 0; i < patch->numAtoms; ++i ) {
      altVirial += outer(patch->f[Results::slow][i],patch->p[i]);
    }
    ADD_TENSOR_OBJECT(reduction,REDUCTION_ALT_VIRIAL_SLOW,altVirial);
  }

  {
    BigReal intKineticEnergy = 0;
    Tensor intVirialNormal;
    Tensor intVirialNbond;
    Tensor intVirialSlow;

    int hgs;
    for ( int i = 0; i < patch->numAtoms; i += hgs ) {
      hgs = patch->a[i].hydrogenGroupSize;
      int j;
      BigReal m_cm = 0;
      Position x_cm(0,0,0);
      Velocity v_cm(0,0,0);
      for ( j = i; j < (i+hgs); ++j ) {
        m_cm += patch->a[j].mass;
        x_cm += patch->a[j].mass * patch->p[j];
        v_cm += patch->a[j].mass * patch->v[j];
      }
      x_cm /= m_cm;
      v_cm /= m_cm;
      for ( j = i; j < (i+hgs); ++j ) {
        BigReal mass = patch->a[j].mass;
        Vector v = patch->v[j];
        Vector dv = v - v_cm;
        intKineticEnergy += mass * (v * dv);
        intVirialNormal += mass * outer(v,dv);
        Vector dx = patch->p[j] - x_cm;
        intVirialNormal += outer(patch->f[Results::normal][j],dx);
        intVirialNbond += outer(patch->f[Results::nbond][j],dx);
        intVirialSlow += outer(patch->f[Results::slow][j],dx);
      }
    }

    intKineticEnergy *= 0.5;

    reduction->item(REDUCTION_INT_KINETIC_ENERGY) += intKineticEnergy;
    ADD_TENSOR_OBJECT(reduction,REDUCTION_INT_VIRIAL_NORMAL,intVirialNormal);
    ADD_TENSOR_OBJECT(reduction,REDUCTION_INT_VIRIAL_NBOND,intVirialNbond);
    ADD_TENSOR_OBJECT(reduction,REDUCTION_INT_VIRIAL_SLOW,intVirialSlow);
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
  int prec = Output::coordinateNeeded(step);
  if ( prec )
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
    Tensor virial;
    patch->mollyMollify(&virial);
    ADD_TENSOR_OBJECT(reduction,REDUCTION_VIRIAL_SLOW,virial);
  }
}

void Sequencer::rebalanceLoad(int timestep) {
  if ( ! ldbSteps ) {
    ldbSteps = LdbCoordinator::Object()->steps();
  }
  if ( ! --ldbSteps ) {
    patch->submitLoadStats(timestep);
    ldbCoordinator->rebalance(this,patch->getPatchID());
  }
}

void Sequencer::cycleBarrier(int doBarrier, int step) {
#ifdef CYCLE_BARRIER
	if (doBarrier)
	  broadcast->cycleBarrier.get(step);
#endif
}

void Sequencer::terminate() {
  CthFree(thread);
  CthSuspend();
}

