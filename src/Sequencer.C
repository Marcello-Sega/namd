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
      case SCRIPT_MINIMIZE:
	minimize();
	continue;
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

    const int commOnly = simParams->commOnly;

    int &maxForceUsed = patch->flags.maxForceUsed;
    int &maxForceMerged = patch->flags.maxForceMerged;
    maxForceUsed = Results::normal;
    maxForceMerged = Results::normal;

    const int numberOfSteps = simParams->N;
    const int stepsPerCycle = simParams->stepsPerCycle;
    const BigReal timestep = simParams->dt;

    // what MTS method?
    const int staleForces = ( simParams->MTSAlgorithm == NAIVE );

    const int nonbondedFrequency = simParams->nonbondedFrequency;
    slowFreq = nonbondedFrequency;
    const BigReal nbondstep = timestep * (staleForces?1:nonbondedFrequency);
    int &doNonbonded = patch->flags.doNonbonded;
    doNonbonded = !(step%nonbondedFrequency);
    if ( nonbondedFrequency == 1 ) maxForceMerged = Results::nbond;
    if ( doNonbonded ) maxForceUsed = Results::nbond;

    // Do we do full electrostatics?
    const int dofull = ( simParams->fullDirectOn ||
			simParams->FMAOn || simParams->PMEOn );
    const int fullElectFrequency = simParams->fullElectFrequency;
    if ( dofull ) slowFreq = fullElectFrequency;
    const BigReal slowstep = timestep * (staleForces?1:fullElectFrequency);
    int &doFullElectrostatics = patch->flags.doFullElectrostatics;
    doFullElectrostatics = (dofull && !(step%fullElectFrequency));
    if ( dofull && (fullElectFrequency == 1) && !(simParams->mollyOn) )
					maxForceMerged = Results::slow;
    if ( doFullElectrostatics ) maxForceUsed = Results::slow;
    int &doMolly = patch->flags.doMolly;
    doMolly = simParams->mollyOn && doFullElectrostatics;

    rattle1(0.);  // enforce rigid bond constraints on initial positions
    minimizationQuenchVelocity();
    runComputeObjects(1); // must migrate here!
    if ( staleForces ) {
      if ( doNonbonded ) saveForce(Results::nbond);
      if ( doFullElectrostatics ) saveForce(Results::slow);
    }
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

       if ( ! commOnly ) {
	addForceToMomentum(0.5*timestep);
	if (staleForces || doNonbonded)
		addForceToMomentum(0.5*nbondstep,Results::nbond,staleForces);
	if (staleForces || doFullElectrostatics)
		addForceToMomentum(0.5*slowstep,Results::slow,staleForces);
	langevinVelocitiesBBK1(timestep);
       }

	maximumMove(timestep);
	if ( ! commOnly ) addVelocityToPosition(0.5*timestep);
	langevinPiston(step);
	if ( ! commOnly ) addVelocityToPosition(0.5*timestep);
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
	if ( staleForces ) {
	  if ( doNonbonded ) saveForce(Results::nbond);
	  if ( doFullElectrostatics ) saveForce(Results::slow);
	}

       if ( ! commOnly ) {
	langevinVelocitiesBBK2(timestep);
	addForceToMomentum(0.5*timestep);
	if (staleForces || doNonbonded)
		addForceToMomentum(0.5*nbondstep,Results::nbond,staleForces);
	if (staleForces || doFullElectrostatics)
		addForceToMomentum(0.5*slowstep,Results::slow,staleForces);
       }

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

void Sequencer::minimize() {
  const int numberOfSteps = simParams->N;
  int &step = patch->flags.step;
  step = simParams->firstTimestep;

  int &maxForceUsed = patch->flags.maxForceUsed;
  int &maxForceMerged = patch->flags.maxForceMerged;
  maxForceUsed = Results::normal;
  maxForceMerged = Results::normal;
  int &doNonbonded = patch->flags.doNonbonded;
  doNonbonded = 1;
  maxForceUsed = Results::nbond;
  maxForceMerged = Results::nbond;
  const int dofull = ( simParams->fullDirectOn ||
			simParams->FMAOn || simParams->PMEOn );
  int &doFullElectrostatics = patch->flags.doFullElectrostatics;
  doFullElectrostatics = dofull;
  if ( dofull ) {
    maxForceMerged = Results::slow;
    maxForceUsed = Results::slow;
  }
  int &doMolly = patch->flags.doMolly;
  doMolly = simParams->mollyOn && doFullElectrostatics;

  rattle1(0.);  // enforce rigid bond constraints on initial positions
  runComputeObjects(1); // must migrate here!

  // addForceToMomentum(0.); // zero velocities of fixed atoms
  submitMinimizeReductions(step);
  rebalanceLoad(step);

  int minSeq = 0;
  for ( ++step; step <= numberOfSteps; ++step ) {
    BigReal c = broadcast->minimizeCoefficient.get(minSeq++);
    if ( ! c ) {  // new direction
      patch->checkpoint();  // checkpoint current minimum
      c = broadcast->minimizeCoefficient.get(minSeq++);
      newMinimizeDirection(c);  // v = c * v + f
      c = broadcast->minimizeCoefficient.get(minSeq++);
    }  // same direction
    newMinimizePosition(c);  // x = x + c * v

    rattle1(0.);
    runComputeObjects(1);
    submitMinimizeReductions(step);
    submitCollections(step);
    rebalanceLoad(step);
  }
  patch->revert();  // go to last known minimum
}

// v = c * v + f
void Sequencer::newMinimizeDirection(BigReal c) {
  FullAtom *a = patch->atom.begin();
  Force *f1 = patch->f[Results::normal].begin();
  Force *f2 = patch->f[Results::nbond].begin();
  Force *f3 = patch->f[Results::slow].begin();
  int numAtoms = patch->numAtoms;

  for ( int i = 0; i < numAtoms; ++i ) {
    a[i].velocity *= c;
    a[i].velocity += f1[i] + f2[i] + f3[i];
  }
}

// x = x + c * v
void Sequencer::newMinimizePosition(BigReal c) {
  FullAtom *a = patch->atom.begin();
  int numAtoms = patch->numAtoms;

  for ( int i = 0; i < numAtoms; ++i ) {
    a[i].position += c * a[i].velocity;
  }
}

void Sequencer::langevinVelocities(BigReal dt_fs)
{
  if ( simParams->langevinOn )
  {
    FullAtom *a = patch->atom.begin();
    int numAtoms = patch->numAtoms;
    Molecule *molecule = Node::Object()->molecule;
    BigReal dt = dt_fs * 0.001;  // convert to ps
    BigReal kbT = BOLTZMAN*(simParams->langevinTemp);
    for ( int i = 0; i < numAtoms; ++i )
    {
      int aid = a[i].id;
      BigReal f1 = exp( -1. * dt * molecule->langevin_param(aid) );
      BigReal f2 = sqrt( ( 1. - f1*f1 ) * kbT / a[i].mass );

      a[i].velocity *= f1;
      a[i].velocity += f2 * random->gaussian_vector();
    }
  }
}

void Sequencer::langevinVelocitiesBBK1(BigReal dt_fs)
{
  if ( simParams->langevinOn )
  {
    FullAtom *a = patch->atom.begin();
    int numAtoms = patch->numAtoms;
    Molecule *molecule = Node::Object()->molecule;
    BigReal dt = dt_fs * 0.001;  // convert to ps
    BigReal kbT = BOLTZMAN*(simParams->langevinTemp);
    for ( int i = 0; i < numAtoms; ++i )
    {
      int aid = a[i].id;
      BigReal dt_gamma = dt * molecule->langevin_param(aid);

      a[i].velocity += random->gaussian_vector() *
             sqrt( dt_gamma * kbT / a[i].mass );
      a[i].velocity /= ( 1. + 0.5 * dt_gamma );
    }
  }
}


void Sequencer::langevinVelocitiesBBK2(BigReal dt_fs)
{
  if ( simParams->langevinOn )
  {
    FullAtom *a = patch->atom.begin();
    int numAtoms = patch->numAtoms;
    Molecule *molecule = Node::Object()->molecule;
    BigReal dt = dt_fs * 0.001;  // convert to ps
    BigReal kbT = BOLTZMAN*(simParams->langevinTemp);
    for ( int i = 0; i < numAtoms; ++i )
    {
      int aid = a[i].id;
      BigReal dt_gamma = dt * molecule->langevin_param(aid);

      a[i].velocity *= ( 1. - 0.5 * dt_gamma );
      a[i].velocity += random->gaussian_vector() *
             sqrt( dt_gamma * kbT / a[i].mass );
    }
  }
}

void Sequencer::berendsenPressure(int step)
{
  const int freq = simParams->berendsenPressureFreq;
  if ( simParams->berendsenPressureOn && !((step-1)%freq) )
  {
   FullAtom *a = patch->atom.begin();
   int numAtoms = patch->numAtoms;
   Tensor factor = broadcast->positionRescaleFactor.get(step);
   patch->lattice.rescale(factor);
   if ( simParams->useGroupPressure )
   {
    int hgs;
    for ( int i = 0; i < numAtoms; i += hgs ) {
      hgs = a[i].hydrogenGroupSize;
      int j;
      BigReal m_cm = 0;
      Position x_cm(0,0,0);
      for ( j = i; j < (i+hgs); ++j ) {
        m_cm += a[j].mass;
        x_cm += a[j].mass * a[j].position;
      }
      x_cm /= m_cm;
      Position new_x_cm = x_cm;
      patch->lattice.rescale(new_x_cm,factor);
      Position delta_x_cm = new_x_cm - x_cm;
      for ( j = i; j < (i+hgs); ++j ) {
        a[j].position += delta_x_cm;
      }
    }
   }
   else
   {
    for ( int i = 0; i < numAtoms; ++i )
    {
      patch->lattice.rescale(a[i].position,factor);
    }
   }
  }
}

void Sequencer::langevinPiston(int step)
{
  if ( simParams->langevinPistonOn && ! ( (step-1-slowFreq/2) % slowFreq ) )
  {
   FullAtom *a = patch->atom.begin();
   int numAtoms = patch->numAtoms;
   Tensor factor = broadcast->positionRescaleFactor.get(step);
   // JCP FIX THIS!!!
   Vector velFactor(1/factor.xx,1/factor.yy,1/factor.zz);
   patch->lattice.rescale(factor);
   if ( simParams->useGroupPressure )
   {
    int hgs;
    for ( int i = 0; i < numAtoms; i += hgs ) {
      hgs = a[i].hydrogenGroupSize;
      int j;
      BigReal m_cm = 0;
      Position x_cm(0,0,0);
      Velocity v_cm(0,0,0);
      for ( j = i; j < (i+hgs); ++j ) {
        m_cm += a[j].mass;
        x_cm += a[j].mass * a[j].position;
        v_cm += a[j].mass * a[j].velocity;
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
        a[j].position += delta_x_cm;
        a[j].velocity += delta_v_cm;
      }
    }
   }
   else
   {
    for ( int i = 0; i < numAtoms; ++i )
    {
      patch->lattice.rescale(a[i].position,factor);
      a[i].velocity.x *= velFactor.x;
      a[i].velocity.y *= velFactor.y;
      a[i].velocity.z *= velFactor.z;
    }
   }
  }
}

void Sequencer::rescaleVelocities(int step)
{
  const int rescaleFreq = simParams->rescaleFreq;
  if ( rescaleFreq > 0 ) {
    FullAtom *a = patch->atom.begin();
    int numAtoms = patch->numAtoms;
    ++rescaleVelocities_numTemps;
    if ( rescaleVelocities_numTemps == rescaleFreq ) {
      BigReal factor = broadcast->velocityRescaleFactor.get(step);
      for ( int i = 0; i < numAtoms; ++i )
      {
        a[i].velocity *= factor;
      }
      rescaleVelocities_numTemps = 0;
    }
  }
}

void Sequencer::reassignVelocities(int step)
{
  const int reassignFreq = simParams->reassignFreq;
  if ( ( reassignFreq > 0 ) && ! ( step % reassignFreq ) ) {
    FullAtom *a = patch->atom.begin();
    int numAtoms = patch->numAtoms;
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
    for ( int i = 0; i < numAtoms; ++i )
    {
      a[i].velocity = ( ( a[i].atomFixed ) ? Vector(0,0,0) :
        sqrt( kbT / a[i].mass ) * random->gaussian_vector() );
    }
  }
}

void Sequencer::reinitVelocities(void)
{
  FullAtom *a = patch->atom.begin();
  int numAtoms = patch->numAtoms;
  BigReal newTemp = simParams->initialTemp;
  BigReal kbT = BOLTZMAN * newTemp;
  for ( int i = 0; i < numAtoms; ++i )
  {
    a[i].velocity = ( ( a[i].atomFixed ) ? Vector(0,0,0) :
      sqrt( kbT / a[i].mass ) * random->gaussian_vector() );
  }
}

void Sequencer::tcoupleVelocities(BigReal dt_fs, int step)
{
  if ( simParams->tCoupleOn )
  {
    FullAtom *a = patch->atom.begin();
    int numAtoms = patch->numAtoms;
    BigReal coefficient = broadcast->tcoupleCoefficient.get(step);
    Molecule *molecule = Node::Object()->molecule;
    BigReal dt = dt_fs * 0.001;  // convert to ps
    coefficient *= dt;
    for ( int i = 0; i < numAtoms; ++i )
    {
      int aid = a[i].id;
      BigReal f1 = exp( coefficient * molecule->langevin_param(aid) );
      a[i].velocity *= f1;
    }
  }
}

void Sequencer::saveForce(const int ftag)
{
  patch->saveForce(ftag);
}

void Sequencer::addForceToMomentum(BigReal dt, const int ftag,
						const int useSaved)
{
  patch->addForceToMomentum(dt,ftag,useSaved);
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
#ifdef ALTVIRIAL
    ADD_TENSOR_OBJECT(reduction,REDUCTION_ALT_VIRIAL_NORMAL,virial);
#endif
    ADD_TENSOR_OBJECT(reduction,REDUCTION_INT_VIRIAL_NORMAL,virial);
  }
}

void Sequencer::maximumMove(BigReal timestep)
{
  FullAtom *a = patch->atom.begin();
  int numAtoms = patch->numAtoms;
  if ( simParams->maximumMove ) {
    const BigReal dt = timestep / TIMEFACTOR;
    const BigReal maxvel = simParams->maximumMove / dt;
    const BigReal maxvel2 = maxvel * maxvel;
    for ( int i=0; i<numAtoms; ++i ) {
      if ( a[i].velocity.length2() > maxvel2 ) {
	a[i].velocity *= ( maxvel / a[i].velocity.length() );
      }
    }
  } else {
    const BigReal dt = timestep / TIMEFACTOR;
    const BigReal maxvel = 10.0 / dt;
    const BigReal maxvel2 = maxvel * maxvel;
    int killme = 0;
    for ( int i=0; i<numAtoms; ++i ) {
      killme = killme || ( a[i].velocity.length2() > maxvel2 );
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
    FullAtom *a = patch->atom.begin();
    int numAtoms = patch->numAtoms;
    for ( int i=0; i<numAtoms; ++i ) {
      a[i].velocity = 0.;
    }
  }
}

void Sequencer::submitReductions(int step)
{
  FullAtom *a = patch->atom.begin();
  int numAtoms = patch->numAtoms;

  reduction->item(REDUCTION_ATOM_CHECKSUM) += numAtoms;
  reduction->item(REDUCTION_KINETIC_ENERGY) += patch->calcKineticEnergy();

  {
    Tensor virial;
    for ( int i = 0; i < numAtoms; ++i ) {
      virial += ( a[i].mass * outer(a[i].velocity,a[i].velocity) );
    }
    ADD_TENSOR_OBJECT(reduction,REDUCTION_VIRIAL_NORMAL,virial);
#ifdef ALTVIRIAL
    ADD_TENSOR_OBJECT(reduction,REDUCTION_ALT_VIRIAL_NORMAL,virial);
  }
  {
    Tensor altVirial;
    for ( int i = 0; i < numAtoms; ++i ) {
      altVirial += outer(patch->f[Results::normal][i],a[i].position);
    }
    ADD_TENSOR_OBJECT(reduction,REDUCTION_ALT_VIRIAL_NORMAL,altVirial);
  }
  {
    Tensor altVirial;
    for ( int i = 0; i < numAtoms; ++i ) {
      altVirial += outer(patch->f[Results::nbond][i],a[i].position);
    }
    ADD_TENSOR_OBJECT(reduction,REDUCTION_ALT_VIRIAL_NBOND,altVirial);
  }
  {
    Tensor altVirial;
    for ( int i = 0; i < numAtoms; ++i ) {
      altVirial += outer(patch->f[Results::slow][i],a[i].position);
    }
    ADD_TENSOR_OBJECT(reduction,REDUCTION_ALT_VIRIAL_SLOW,altVirial);
#endif
  }

  {
    BigReal intKineticEnergy = 0;
    Tensor intVirialNormal;
    Tensor intVirialNbond;
    Tensor intVirialSlow;

    int hgs;
    for ( int i = 0; i < numAtoms; i += hgs ) {
      hgs = a[i].hydrogenGroupSize;
      int j;
      BigReal m_cm = 0;
      Position x_cm(0,0,0);
      Velocity v_cm(0,0,0);
      for ( j = i; j < (i+hgs); ++j ) {
        m_cm += a[j].mass;
        x_cm += a[j].mass * a[j].position;
        v_cm += a[j].mass * a[j].velocity;
      }
      x_cm /= m_cm;
      v_cm /= m_cm;
      for ( j = i; j < (i+hgs); ++j ) {
        BigReal mass = a[j].mass;
        Vector v = a[j].velocity;
        Vector dv = v - v_cm;
        intKineticEnergy += mass * (v * dv);
        intVirialNormal += mass * outer(v,dv);
        Vector dx = a[j].position - x_cm;
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

void Sequencer::submitMinimizeReductions(int step)
{
  FullAtom *a = patch->atom.begin();
  Force *f1 = patch->f[Results::normal].begin();
  Force *f2 = patch->f[Results::nbond].begin();
  Force *f3 = patch->f[Results::slow].begin();
  int numAtoms = patch->numAtoms;

  reduction->item(REDUCTION_ATOM_CHECKSUM) += numAtoms;
  const BigReal fmax = 100. * TIMEFACTOR * TIMEFACTOR;
  const BigReal fmax2 = fmax * fmax;

  BigReal fdotf = 0;
  BigReal fdotv = 0;
  BigReal vdotv = 0;
  for ( int i = 0; i < numAtoms; ++i ) {
    Force f = f1[i] + f2[i] + f3[i];
    BigReal ff = f * f;
    if ( ff > fmax2 ) {
      BigReal fmult = sqrt(fmax2/ff);
      f *= fmult;  ff = f * f;
      f1[i] *= fmult;
      f2[i] *= fmult;
      f3[i] *= fmult;
    }
    fdotf += ff;
    fdotv += f * a[i].velocity;
    vdotv += a[i].velocity * a[i].velocity;
  }

  reduction->item(REDUCTION_MIN_F_DOT_F) += fdotf;
  reduction->item(REDUCTION_MIN_F_DOT_V) += fdotv;
  reduction->item(REDUCTION_MIN_V_DOT_V) += vdotv;

  reduction->submit();
}

void Sequencer::submitCollections(int step)
{
  int prec = Output::coordinateNeeded(step);
  if ( prec )
    collection->submitPositions(step,patch->atom,patch->lattice,prec);
  if ( Output::velocityNeeded(step) )
    collection->submitVelocities(step,patch->atom);
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

