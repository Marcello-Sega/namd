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
    berendsenPressure_count = 0;
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
  int scriptTask;
  int scriptSeq = 0;
  while ( (scriptTask = broadcast->scriptBarrier.get(scriptSeq++)) != SCRIPT_END ) {
    switch ( scriptTask ) {
      case SCRIPT_OUTPUT:
	submitCollections(FILE_OUTPUT);
	break;
      case SCRIPT_MEASURE:
	submitCollections(EVAL_MEASURE);
	break;
      case SCRIPT_REINITVELS:
	reinitVelocities();
	break;
      case SCRIPT_CHECKPOINT:
        patch->checkpoint();
        checkpoint_berendsenPressure_count = berendsenPressure_count;
	break;
      case SCRIPT_REVERT:
        patch->revert();
        berendsenPressure_count = checkpoint_berendsenPressure_count;
	break;
      case SCRIPT_MINIMIZE:
	minimize();
	break;
      case SCRIPT_RUN:
	integrate();
	break;
    }
  }
  submitCollections(END_OF_RUN);
  terminate();
}


void Sequencer::integrate() {

    int &step = patch->flags.step;
    step = simParams->firstTimestep;

    // drag switches
    const Bool rotDragOn = simParams->rotDragOn;
    const Bool movDragOn = simParams->movDragOn;

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

        // add drag to each atom's positions
        if ( ! commOnly && movDragOn ) addMovDragToPosition(timestep);
        if ( ! commOnly && rotDragOn ) addRotDragToPosition(timestep);

	submitReductions(step);
	submitCollections(step);
        cycleBarrier(!((step+1) % stepsPerCycle),step);
	rebalanceLoad(step);
    }
}

// add moving drag to each atom's position
void Sequencer::addMovDragToPosition(BigReal timestep) {
  FullAtom *atom = patch->atom.begin();
  int numAtoms = patch->numAtoms;
  Molecule *molecule = Node::Object()->molecule;   // need its methods
  const BigReal movDragGlobVel = simParams->movDragGlobVel;
  const BigReal dt = timestep / TIMEFACTOR;   // MUST be as in the integrator!
  Vector movDragVel, dragIncrement;
  for ( int i = 0; i < numAtoms; ++i )
  {
    // skip if fixed atom or zero drag attribute
    if ( (simParams->fixedAtomsOn && atom[i].atomFixed) 
	 || !(molecule->is_atom_movdragged(atom[i].id)) ) continue;
    molecule->get_movdrag_params(movDragVel, atom[i].id);
    dragIncrement = movDragGlobVel * movDragVel * dt;
    atom[i].position += dragIncrement;
  }
}

// add rotating drag to each atom's position
void Sequencer::addRotDragToPosition(BigReal timestep) {
  FullAtom *atom = patch->atom.begin();
  int numAtoms = patch->numAtoms;
  Molecule *molecule = Node::Object()->molecule;   // need its methods
  const BigReal rotDragGlobVel = simParams->rotDragGlobVel;
  const BigReal dt = timestep / TIMEFACTOR;   // MUST be as in the integrator!
  BigReal rotDragVel, dAngle;
  Vector atomRadius, atomTangent, atomNormal;
  Vector rotDragUnit;
  Vector rotDragAxis, rotDragPivot, dragIncrement;
  for ( int i = 0; i < numAtoms; ++i )
  {
    // skip if fixed atom or zero drag attribute
    if ( (simParams->fixedAtomsOn && atom[i].atomFixed) 
	 || !(molecule->is_atom_rotdragged(atom[i].id)) ) continue;
    molecule->get_rotdrag_params(rotDragVel, rotDragAxis, rotDragPivot, atom[i].id);
    dAngle = rotDragGlobVel * rotDragVel * dt;
    rotDragUnit = rotDragAxis / rotDragAxis.length();
    atomRadius = atom[i].position - rotDragPivot;
    atomTangent = rotDragUnit * (atomRadius * rotDragUnit) / atomRadius.length();
    atomNormal = atomRadius - atomTangent;
    dragIncrement = cross(rotDragUnit, atomNormal) * dAngle;
    atom[i].position += dragIncrement;
  }
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

  runComputeObjects(1); // must migrate here!

  submitMinimizeReductions(step);
  rebalanceLoad(step);

  int minSeq = 0;
  for ( ++step; step <= numberOfSteps; ++step ) {
    BigReal c = broadcast->minimizeCoefficient.get(minSeq++);
    if ( ! c ) {  // new direction
      c = broadcast->minimizeCoefficient.get(minSeq++);
      newMinimizeDirection(c);  // v = c * v + f
      c = broadcast->minimizeCoefficient.get(minSeq++);
    }  // same direction
    newMinimizePosition(c);  // x = x + c * v

    runComputeObjects(1);
    submitMinimizeReductions(step);
    submitCollections(step);
    rebalanceLoad(step);
  }
  quenchVelocities();  // zero out bogus velocity
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
    if ( simParams->fixedAtomsOn && a[i].atomFixed ) a[i].velocity = 0;
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

// v = 0
void Sequencer::quenchVelocities() {
  FullAtom *a = patch->atom.begin();
  int numAtoms = patch->numAtoms;

  for ( int i = 0; i < numAtoms; ++i ) {
    a[i].velocity = 0;
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
      if ( ! dt_gamma ) continue;

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
      if ( ! dt_gamma ) continue;

      a[i].velocity *= ( 1. - 0.5 * dt_gamma );
      a[i].velocity += random->gaussian_vector() *
             sqrt( dt_gamma * kbT / a[i].mass );
    }
  }
}

void Sequencer::berendsenPressure(int step)
{
  if ( simParams->berendsenPressureOn ) {
  berendsenPressure_count += 1;
  const int freq = simParams->berendsenPressureFreq;
  if ( ! (berendsenPressure_count % freq ) ) {
   berendsenPressure_count = 0;
   FullAtom *a = patch->atom.begin();
   int numAtoms = patch->numAtoms;
   Tensor factor = broadcast->positionRescaleFactor.get(step);
   patch->lattice.rescale(factor);
   if ( simParams->useGroupPressure )
   {
    int hgs;
    for ( int i = 0; i < numAtoms; i += hgs ) {
      int j;
      hgs = a[i].hydrogenGroupSize;
      if ( simParams->fixedAtomsOn && a[i].groupFixed ) {
        for ( j = i; j < (i+hgs); ++j ) {
          a[j].position = patch->lattice.apply_transform(
				a[j].fixedPosition,a[j].transform);
        }
        continue;
      }
      BigReal m_cm = 0;
      Position x_cm(0,0,0);
      for ( j = i; j < (i+hgs); ++j ) {
        if ( simParams->fixedAtomsOn && a[j].atomFixed ) continue;
        m_cm += a[j].mass;
        x_cm += a[j].mass * a[j].position;
      }
      x_cm /= m_cm;
      Position new_x_cm = x_cm;
      patch->lattice.rescale(new_x_cm,factor);
      Position delta_x_cm = new_x_cm - x_cm;
      for ( j = i; j < (i+hgs); ++j ) {
        if ( simParams->fixedAtomsOn && a[j].atomFixed ) {
          a[j].position = patch->lattice.apply_transform(
				a[j].fixedPosition,a[j].transform);
          continue;
        }
        a[j].position += delta_x_cm;
      }
    }
   }
   else
   {
    for ( int i = 0; i < numAtoms; ++i )
    {
      if ( simParams->fixedAtomsOn && a[i].atomFixed ) {
        a[i].position = patch->lattice.apply_transform(
				a[i].fixedPosition,a[i].transform);
        continue;
      }
      patch->lattice.rescale(a[i].position,factor);
    }
   }
  }
  } else {
    berendsenPressure_count = 0;
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
   Molecule *mol = Node::Object()->molecule;
   if ( simParams->useGroupPressure )
   {
    int hgs;
    for ( int i = 0; i < numAtoms; i += hgs ) {
      int j;
      hgs = a[i].hydrogenGroupSize;
      if ( simParams->fixedAtomsOn && a[i].groupFixed ) {
        for ( j = i; j < (i+hgs); ++j ) {
          a[j].position = patch->lattice.apply_transform(
				a[j].fixedPosition,a[j].transform);
        }
        continue;
      }
      BigReal m_cm = 0;
      Position x_cm(0,0,0);
      Velocity v_cm(0,0,0);
      for ( j = i; j < (i+hgs); ++j ) {
        if ( simParams->fixedAtomsOn && a[j].atomFixed ) continue;
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
        if ( simParams->fixedAtomsOn && a[j].atomFixed ) {
          a[j].position = patch->lattice.apply_transform(
				a[j].fixedPosition,a[j].transform);
          continue;
        }
        if ( mol->is_atom_exPressure(a[j].id) ) continue;
        a[j].position += delta_x_cm;
        a[j].velocity += delta_v_cm;
      }
    }
   }
   else
   {
    for ( int i = 0; i < numAtoms; ++i )
    {
      if ( simParams->fixedAtomsOn && a[i].atomFixed ) {
        a[i].position = patch->lattice.apply_transform(
				a[i].fixedPosition,a[i].transform);
        continue;
      }
      if ( mol->is_atom_exPressure(a[i].id) ) continue;
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
      a[i].velocity = ( ( simParams->fixedAtomsOn && a[i].atomFixed ) ? Vector(0,0,0) :
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
    a[i].velocity = ( ( simParams->fixedAtomsOn && a[i].atomFixed ) ? Vector(0,0,0) :
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
    if ( patch->rattle1(dt) ) {
      iout << iERROR << 
        "Constraint failure; simulation has become unstable.\n" << endi;
      Node::Object()->enableEarlyExit();
      terminate();
    }
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
      for ( int i=0; i<numAtoms; ++i ) {
        if ( a[i].velocity.length2() > maxvel2 ) {
          iout << iERROR << "Atom " << (a[i].id + 1) << " velocity is "
            << ( PDBVELFACTOR * a[i].velocity ) << " (limit is "
            << ( PDBVELFACTOR * maxvel ) << ")\n" << endi;
        }
      }
      iout << iERROR << 
        "Atoms moving too fast; simulation has become unstable.\n" << endi;
      Node::Object()->enableEarlyExit();
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
  reduction->item(REDUCTION_MARGIN_VIOLATIONS) += patch->marginViolations;
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
	// net force treated as zero for fixed atoms
        if ( simParams->fixedAtomsOn && a[j].atomFixed ) continue;
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

  if ( simParams->fixedAtomsOn ) {
    Tensor fixVirialNormal;
    Tensor fixVirialNbond;
    Tensor fixVirialSlow;
    Vector fixForceNormal = 0;
    Vector fixForceNbond = 0;
    Vector fixForceSlow = 0;
    Vector fixPosition = 0;
    int fixCount = 0;

    for ( int j = 0; j < numAtoms; j++ ) {
      if ( simParams->fixedAtomsOn && a[j].atomFixed ) {
        Vector dx = a[j].fixedPosition;
        // all negative because fixed atoms cancels these forces
        fixVirialNormal -= outer(patch->f[Results::normal][j],dx);
        fixVirialNbond -= outer(patch->f[Results::nbond][j],dx);
        fixVirialSlow -= outer(patch->f[Results::slow][j],dx);
        fixForceNormal -= patch->f[Results::normal][j];
        fixForceNbond -= patch->f[Results::nbond][j];
        fixForceSlow -= patch->f[Results::slow][j];
        fixPosition += dx;
        fixCount += 1;
      }
    }

    ADD_TENSOR_OBJECT(reduction,REDUCTION_VIRIAL_NORMAL,fixVirialNormal);
    ADD_TENSOR_OBJECT(reduction,REDUCTION_VIRIAL_NBOND,fixVirialNbond);
    ADD_TENSOR_OBJECT(reduction,REDUCTION_VIRIAL_SLOW,fixVirialSlow);
    ADD_VECTOR_OBJECT(reduction,REDUCTION_EXT_FORCE_NORMAL,fixForceNormal);
    ADD_VECTOR_OBJECT(reduction,REDUCTION_EXT_FORCE_NBOND,fixForceNbond);
    ADD_VECTOR_OBJECT(reduction,REDUCTION_EXT_FORCE_SLOW,fixForceSlow);
    ADD_VECTOR_OBJECT(reduction,REDUCTION_EXT_POSITION,fixPosition);
    reduction->item(REDUCTION_EXT_COUNT) += fixCount;
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
  int numHuge = 0;
  for ( int i = 0; i < numAtoms; ++i ) {
    if ( simParams->fixedAtomsOn && a[i].atomFixed ) continue;
    Force f = f1[i] + f2[i] + f3[i];
    BigReal ff = f * f;
    if ( ff > fmax2 ) {
      ++numHuge;
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
  reduction->item(REDUCTION_MIN_HUGE_COUNT) += numHuge;

  {
    Tensor intVirialNormal;
    Tensor intVirialNbond;
    Tensor intVirialSlow;

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
      for ( j = i; j < (i+hgs); ++j ) {
        BigReal mass = a[j].mass;
	// net force treated as zero for fixed atoms
        if ( simParams->fixedAtomsOn && a[j].atomFixed ) continue;
        Vector dx = a[j].position - x_cm;
        intVirialNormal += outer(patch->f[Results::normal][j],dx);
        intVirialNbond += outer(patch->f[Results::nbond][j],dx);
        intVirialSlow += outer(patch->f[Results::slow][j],dx);
      }
    }

    ADD_TENSOR_OBJECT(reduction,REDUCTION_INT_VIRIAL_NORMAL,intVirialNormal);
    ADD_TENSOR_OBJECT(reduction,REDUCTION_INT_VIRIAL_NBOND,intVirialNbond);
    ADD_TENSOR_OBJECT(reduction,REDUCTION_INT_VIRIAL_SLOW,intVirialSlow);
  }

  if ( simParams->fixedAtomsOn ) {
    Tensor fixVirialNormal;
    Tensor fixVirialNbond;
    Tensor fixVirialSlow;
    Vector fixForceNormal = 0;
    Vector fixForceNbond = 0;
    Vector fixForceSlow = 0;
    Vector fixPosition = 0;
    int fixCount = 0;

    for ( int j = 0; j < numAtoms; j++ ) {
      if ( simParams->fixedAtomsOn && a[j].atomFixed ) {
        Vector dx = a[j].fixedPosition;
        // all negative because fixed atoms cancels these forces
        fixVirialNormal -= outer(patch->f[Results::normal][j],dx);
        fixVirialNbond -= outer(patch->f[Results::nbond][j],dx);
        fixVirialSlow -= outer(patch->f[Results::slow][j],dx);
        fixForceNormal -= patch->f[Results::normal][j];
        fixForceNbond -= patch->f[Results::nbond][j];
        fixForceSlow -= patch->f[Results::slow][j];
        fixPosition += dx;
        fixCount += 1;
      }
    }

    ADD_TENSOR_OBJECT(reduction,REDUCTION_VIRIAL_NORMAL,fixVirialNormal);
    ADD_TENSOR_OBJECT(reduction,REDUCTION_VIRIAL_NBOND,fixVirialNbond);
    ADD_TENSOR_OBJECT(reduction,REDUCTION_VIRIAL_SLOW,fixVirialSlow);
    ADD_VECTOR_OBJECT(reduction,REDUCTION_EXT_FORCE_NORMAL,fixForceNormal);
    ADD_VECTOR_OBJECT(reduction,REDUCTION_EXT_FORCE_NBOND,fixForceNbond);
    ADD_VECTOR_OBJECT(reduction,REDUCTION_EXT_FORCE_SLOW,fixForceSlow);
    ADD_VECTOR_OBJECT(reduction,REDUCTION_EXT_POSITION,fixPosition);
    reduction->item(REDUCTION_EXT_COUNT) += fixCount;
  }

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

