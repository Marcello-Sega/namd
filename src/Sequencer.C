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
    if (simParams->pressureProfileOn && !simParams->pressureProfileNonbonded) {
      pressureProfileReduction = 
        ReductionMgr::Object()->willSubmit(REDUCTIONS_USER1);
    } else {
      pressureProfileReduction = NULL;
    }
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
    delete pressureProfileReduction;
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
        pairlistsAreValid = 0;
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

    // Bother to calculate energies?
    int &doEnergy = patch->flags.doEnergy;
    int energyFrequency = simParams->outputEnergies;

    rattle1(0.,0);  // enforce rigid bond constraints on initial positions
    doEnergy = ! ( step % energyFrequency );
    runComputeObjects(1,step<numberOfSteps); // must migrate here!
    if ( staleForces ) {
      if ( doNonbonded ) saveForce(Results::nbond);
      if ( doFullElectrostatics ) saveForce(Results::slow);
    }
    if ( ! commOnly ) {
      addForceToMomentum(-0.5*timestep);
      if (staleForces || doNonbonded)
		addForceToMomentum(-0.5*nbondstep,Results::nbond,staleForces);
      if (staleForces || doFullElectrostatics)
		addForceToMomentum(-0.5*slowstep,Results::slow,staleForces);
    }
    minimizationQuenchVelocity();
    rattle1(-timestep,0);
    submitHalfstep(step);
    if ( ! commOnly ) {
      addForceToMomentum(timestep);
      if (staleForces || doNonbonded)
		addForceToMomentum(nbondstep,Results::nbond,staleForces);
      if (staleForces || doFullElectrostatics)
		addForceToMomentum(slowstep,Results::slow,staleForces);
    }
    rattle1(timestep,1);
    submitHalfstep(step);
    if ( ! commOnly ) {
      addForceToMomentum(-0.5*timestep);
      if (staleForces || doNonbonded)
		addForceToMomentum(-0.5*nbondstep,Results::nbond,staleForces);
      if (staleForces || doFullElectrostatics)
		addForceToMomentum(-0.5*slowstep,Results::slow,staleForces);
    }
    submitReductions(step);
    rebalanceLoad(step);

    for ( ++step; step <= numberOfSteps; ++step )
    {
	rescaleVelocities(step);
	tcoupleVelocities(timestep,step);
	berendsenPressure(step);

       if ( ! commOnly ) {
	addForceToMomentum(0.5*timestep);
	if (staleForces || doNonbonded)
		addForceToMomentum(0.5*nbondstep,Results::nbond,staleForces);
	if (staleForces || doFullElectrostatics)
		addForceToMomentum(0.5*slowstep,Results::slow,staleForces);
       }

       const int reassignFreq = simParams->reassignFreq;
       if ( !commOnly && ( reassignFreq>0 ) && ! (step%reassignFreq) ) {
	addVelocityToPosition(0.5*timestep);
        reassignVelocities(timestep,step);
	addVelocityToPosition(0.5*timestep);
	rattle1(0.,0);
	rattle1(-timestep,0);
	addVelocityToPosition(-1.0*timestep);
	rattle1(timestep,0);
       }

	maximumMove(timestep);
	if ( ! commOnly ) addVelocityToPosition(0.5*timestep);
	langevinPiston(step);
	if ( ! commOnly ) addVelocityToPosition(0.5*timestep);

	minimizationQuenchVelocity();
	submitHalfstep(step);

	doNonbonded = !(step%nonbondedFrequency);
	doFullElectrostatics = (dofull && !(step%fullElectFrequency));
	doMolly = simParams->mollyOn && doFullElectrostatics;

        maxForceUsed = Results::normal;
	if ( doNonbonded ) maxForceUsed = Results::nbond;
	if ( doFullElectrostatics ) maxForceUsed = Results::slow;

	// Migrate Atoms on stepsPerCycle
        doEnergy = ! ( step % energyFrequency );
	runComputeObjects(!(step%stepsPerCycle),step<numberOfSteps);
	if ( staleForces ) {
	  if ( doNonbonded ) saveForce(Results::nbond);
	  if ( doFullElectrostatics ) saveForce(Results::slow);
	}

       if ( ! commOnly ) {
	langevinVelocitiesBBK2(timestep);
	addForceToMomentum(timestep);
	if (staleForces || doNonbonded)
		addForceToMomentum(nbondstep,Results::nbond,staleForces);
	if (staleForces || doFullElectrostatics)
		addForceToMomentum(slowstep,Results::slow,staleForces);
	langevinVelocitiesBBK1(timestep);
       }

        // add drag to each atom's positions
        if ( ! commOnly && movDragOn ) addMovDragToPosition(timestep);
        if ( ! commOnly && rotDragOn ) addRotDragToPosition(timestep);

	rattle1(timestep,1);

	submitHalfstep(step);

       if ( ! commOnly ) {
	addForceToMomentum(-0.5*timestep);
	if (staleForces || doNonbonded)
		addForceToMomentum(-0.5*nbondstep,Results::nbond,staleForces);
	if (staleForces || doFullElectrostatics)
		addForceToMomentum(-0.5*slowstep,Results::slow,staleForces);
       }

	// rattle2(timestep,step);

	submitReductions(step);
	submitCollections(step);

#if CYCLE_BARRIER
        cycleBarrier(!((step+1) % stepsPerCycle), step);
#elif PME_BARRIER
        cycleBarrier(doFullElectrostatics, step);
#endif

	rebalanceLoad(step);

#if PME_BARRIER
	// a step before PME
        cycleBarrier(dofull && !((step+1)%fullElectFrequency),step);
#endif
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
  Vector atomRadius;
  Vector rotDragAxis, rotDragPivot, dragIncrement;
  for ( int i = 0; i < numAtoms; ++i )
  {
    // skip if fixed atom or zero drag attribute
    if ( (simParams->fixedAtomsOn && atom[i].atomFixed) 
	 || !(molecule->is_atom_rotdragged(atom[i].id)) ) continue;
    molecule->get_rotdrag_params(rotDragVel, rotDragAxis, rotDragPivot, atom[i].id);
    dAngle = rotDragGlobVel * rotDragVel * dt;
    rotDragAxis /= rotDragAxis.length();
    atomRadius = atom[i].position - rotDragPivot;
    dragIncrement = cross(rotDragAxis, atomRadius) * dAngle;
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
  int &doEnergy = patch->flags.doEnergy;
  doEnergy = 1;

  runComputeObjects(1,0); // must migrate here!

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

    runComputeObjects(1,0);
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
             sqrt( 2 * dt_gamma * kbT / a[i].mass );
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
    // BigReal kbT = BOLTZMAN*(simParams->langevinTemp);
    for ( int i = 0; i < numAtoms; ++i )
    {
      int aid = a[i].id;
      BigReal dt_gamma = dt * molecule->langevin_param(aid);
      if ( ! dt_gamma ) continue;

      a[i].velocity *= ( 1. - 0.5 * dt_gamma );
      // a[i].velocity += random->gaussian_vector() *
      //        sqrt( dt_gamma * kbT / a[i].mass );
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

void Sequencer::reassignVelocities(BigReal timestep, int step)
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
  } else {
    NAMD_bug("Sequencer::reassignVelocities called improperly!");
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

void Sequencer::rattle1(BigReal dt, int pressure)
{
  if ( simParams->rigidBonds != RIGID_NONE ) {
    Tensor virial;
    Tensor *vp = ( pressure ? &virial : 0 );
    if ( patch->rattle1(dt, vp, pressureProfileReduction) ) {
      iout << iERROR << 
        "Constraint failure; simulation has become unstable.\n" << endi;
      Node::Object()->enableEarlyExit();
      terminate();
    }
    ADD_TENSOR_OBJECT(reduction,REDUCTION_VIRIAL_NORMAL,virial);
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

void Sequencer::submitHalfstep(int step)
{
  // velocity-dependent quantities *** ONLY ***
  // positions are not at half-step when called
  FullAtom *a = patch->atom.begin();
  int numAtoms = patch->numAtoms;

  {
    BigReal kineticEnergy = 0;
    Tensor virial;
    for ( int i = 0; i < numAtoms; ++i ) {
      kineticEnergy += a[i].mass * a[i].velocity.length2();
      virial += ( a[i].mass * outer(a[i].velocity,a[i].velocity) );
    }

    kineticEnergy *= 0.5 * 0.5;
    reduction->item(REDUCTION_KINETIC_ENERGY) += kineticEnergy;
    virial *= 0.5;
    ADD_TENSOR_OBJECT(reduction,REDUCTION_VIRIAL_NORMAL,virial);
#ifdef ALTVIRIAL
    ADD_TENSOR_OBJECT(reduction,REDUCTION_ALT_VIRIAL_NORMAL,virial);
#endif
  }
 
  if (pressureProfileReduction) {
    BigReal idz = 1.0/simParams->pressureProfileThickness;
    BigReal zmin = simParams->pressureProfileMin;
    int nslabs = simParams->pressureProfileSlabs;
    int useGroupPressure = simParams->useGroupPressure;

    // Compute kinetic energy partition, possibly subtracting off
    // internal kinetic energy if group pressure is enabled.
    // Since the regular pressure is 1/2 mvv and the internal kinetic
    // term that is subtracted off for the group pressure is
    // 1/2 mv (v-v_cm), the group pressure kinetic contribution is
    // 1/2 m * v * v_cm.  The factor of 1/2 is because submitHalfstep
    // gets called twice per timestep.
    int hgs;
    for (int i=0; i<numAtoms; i += hgs) {
      int j, slab;
      hgs = a[i].hydrogenGroupSize;
      BigReal m_cm = 0;
      Velocity v_cm(0,0,0);
      for (j=i; j< i+hgs; ++j) {
        m_cm += a[j].mass;
        v_cm += a[j].mass * a[j].velocity;
      }
      v_cm /= m_cm;
      for (j=i; j < i+hgs; ++j) {
        BigReal mass = a[j].mass;
        if (! (useGroupPressure && j != i)) {
          BigReal z = a[j].position.z;
          slab = (int)floor((z-zmin)*idz);
          if (slab < 0) slab += nslabs;
          else if (slab >= nslabs) slab -= nslabs;
        }
        BigReal wxx, wyy, wzz;
        if (useGroupPressure) {
          wxx = 0.5*mass * a[j].velocity.x * v_cm.x;
          wyy = 0.5*mass * a[j].velocity.y * v_cm.y;
          wzz = 0.5*mass * a[j].velocity.z * v_cm.z;
        } else {
          wxx = 0.5*mass * a[j].velocity.x * a[j].velocity.x;
          wyy = 0.5*mass * a[j].velocity.y * a[j].velocity.y;
          wzz = 0.5*mass * a[j].velocity.z * a[j].velocity.z;
        }
        pressureProfileReduction->item(3*slab  ) += wxx;
        pressureProfileReduction->item(3*slab+1) += wyy;
        pressureProfileReduction->item(3*slab+2) += wzz;
      }
    }
  } 

  {
    Tensor intVirialNormal;

    int hgs;
    for ( int i = 0; i < numAtoms; i += hgs ) {
      hgs = a[i].hydrogenGroupSize;
      int j;
      BigReal m_cm = 0;
      Velocity v_cm(0,0,0);
      for ( j = i; j < (i+hgs); ++j ) {
        m_cm += a[j].mass;
        v_cm += a[j].mass * a[j].velocity;
      }
      v_cm /= m_cm;
      for ( j = i; j < (i+hgs); ++j ) {
        BigReal mass = a[j].mass;
        Vector v = a[j].velocity;
        Vector dv = v - v_cm;
        intVirialNormal += mass * outer(v,dv);
      }
    }

    intVirialNormal *= 0.5;
    ADD_TENSOR_OBJECT(reduction,REDUCTION_INT_VIRIAL_NORMAL,intVirialNormal);
  }

}

void Sequencer::submitReductions(int step)
{
  FullAtom *a = patch->atom.begin();
  int numAtoms = patch->numAtoms;

  reduction->item(REDUCTION_ATOM_CHECKSUM) += numAtoms;
  reduction->item(REDUCTION_MARGIN_VIOLATIONS) += patch->marginViolations;
  reduction->item(REDUCTION_CENTERED_KINETIC_ENERGY) += patch->calcKineticEnergy();

#ifdef ALTVIRIAL
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
  }
#endif

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
	// net force treated as zero for fixed atoms
        if ( simParams->fixedAtomsOn && a[j].atomFixed ) continue;
        BigReal mass = a[j].mass;
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

  if (pressureProfileReduction && simParams->useGroupPressure) {
    // subtract off internal virial term, calculated as for intVirial.
    BigReal idz = 1.0/simParams->pressureProfileThickness;
    BigReal zmin = simParams->pressureProfileMin;
    int nslabs = simParams->pressureProfileSlabs;
    int useGroupPressure = simParams->useGroupPressure;

    int hgs;
    for (int i=0; i<numAtoms; i += hgs) {
      int j;
      hgs = a[i].hydrogenGroupSize;
      BigReal m_cm = 0;
      Position x_cm(0,0,0);
      for (j=i; j< i+hgs; ++j) {
        m_cm += a[j].mass;
        x_cm += a[j].mass * a[j].position;
      }
      x_cm /= m_cm;
      
      BigReal z = a[i].position.z;
      int slab = (int)floor((z-zmin)*idz);
      if (slab < 0) slab += nslabs;
      else if (slab >= nslabs) slab -= nslabs;
      for (j=i; j < i+hgs; ++j) {
        BigReal mass = a[j].mass;
        Vector dx = a[j].position - x_cm;
        const Vector &fnormal = patch->f[Results::normal][j];
        const Vector &fnbond  = patch->f[Results::nbond][j];
        const Vector &fslow   = patch->f[Results::slow][j];
        BigReal wxx = (fnormal.x + fnbond.x + fslow.x) * dx.x;
        BigReal wyy = (fnormal.y + fnbond.y + fslow.y) * dx.y;
        BigReal wzz = (fnormal.z + fnbond.z + fslow.z) * dx.z;
        pressureProfileReduction->item(3*slab  ) -= wxx;
        pressureProfileReduction->item(3*slab+1) -= wyy;
        pressureProfileReduction->item(3*slab+2) -= wzz;
      }
    }
  }

  if ( simParams->fixedAtomsOn ) {
    Tensor fixVirialNormal;
    Tensor fixVirialNbond;
    Tensor fixVirialSlow;
    Vector fixForceNormal = 0;
    Vector fixForceNbond = 0;
    Vector fixForceSlow = 0;

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
      }
    }

    ADD_TENSOR_OBJECT(reduction,REDUCTION_VIRIAL_NORMAL,fixVirialNormal);
    ADD_TENSOR_OBJECT(reduction,REDUCTION_VIRIAL_NBOND,fixVirialNbond);
    ADD_TENSOR_OBJECT(reduction,REDUCTION_VIRIAL_SLOW,fixVirialSlow);
    ADD_VECTOR_OBJECT(reduction,REDUCTION_EXT_FORCE_NORMAL,fixForceNormal);
    ADD_VECTOR_OBJECT(reduction,REDUCTION_EXT_FORCE_NBOND,fixForceNbond);
    ADD_VECTOR_OBJECT(reduction,REDUCTION_EXT_FORCE_SLOW,fixForceSlow);
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
  if (pressureProfileReduction) pressureProfileReduction->submit();
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
      }
    }

    ADD_TENSOR_OBJECT(reduction,REDUCTION_VIRIAL_NORMAL,fixVirialNormal);
    ADD_TENSOR_OBJECT(reduction,REDUCTION_VIRIAL_NBOND,fixVirialNbond);
    ADD_TENSOR_OBJECT(reduction,REDUCTION_VIRIAL_SLOW,fixVirialSlow);
    ADD_VECTOR_OBJECT(reduction,REDUCTION_EXT_FORCE_NORMAL,fixForceNormal);
    ADD_VECTOR_OBJECT(reduction,REDUCTION_EXT_FORCE_NBOND,fixForceNbond);
    ADD_VECTOR_OBJECT(reduction,REDUCTION_EXT_FORCE_SLOW,fixForceSlow);
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

void Sequencer::runComputeObjects(int migration, int pairlists)
{
  if ( migration ) pairlistsAreValid = 0;
  if ( ! simParams->usePairlists ) pairlists = 0;
  patch->flags.usePairlists = pairlists || pairlistsAreValid;
  patch->flags.savePairlists =
	pairlists && ! pairlistsAreValid;
  patch->positionsReady(migration);
  suspend(); // until all deposit boxes close
  if ( patch->flags.savePairlists ) pairlistsAreValid = 1;
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
    pairlistsAreValid = 0;
  }
}

void Sequencer::cycleBarrier(int doBarrier, int step) {
#if USE_BARRIER
	if (doBarrier)
	  broadcast->cycleBarrier.get(step);
#endif
}

void Sequencer::terminate() {
  CthFree(thread);
  CthSuspend();
}

