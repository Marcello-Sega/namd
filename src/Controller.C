/**
***  Copyright (c) 1995, 1996, 1997, 1998, 1999, 2000 by
***  The Board of Trustees of the University of Illinois.
***  All rights reserved.
**/

#include "memusage.h"
#include "Node.h"
#include "Molecule.h"
#include "SimParameters.h"
#include "Controller.h"
#include "ReductionMgr.h"
#include "CollectionMaster.h"
#include "Output.h"
#include "strlib.h"
#include "BroadcastObject.h"
#include "NamdState.h"
#include "ScriptTcl.h"
#include "Broadcasts.h"
#include "LdbCoordinator.h"
#include "Thread.h"
#include <math.h>
#include "NamdOneTools.h"
#include "PatchMap.h"
#include "PatchMap.inl"
#include "Random.h"
#include "imd.h"
#include "IMDOutput.h"
#include "InfoStream.h"

#if(CMK_CCS_AVAILABLE)
extern "C" void CApplicationDepositNode0Data(char *);
#endif

#ifndef cbrt
  // cbrt() not in math.h on goneril
  #define cbrt(x)  pow(x,(double)(1.0/3.0))
#endif

//#define DEBUG_PRESSURE
#define MIN_DEBUG_LEVEL 3
//#define DEBUGM
#include "Debug.h"

Controller::Controller(NamdState *s) :
	computeChecksum(0),
	simParams(Node::Object()->simParameters),
	state(s),
	scriptSeq(0),
	collection(CollectionMaster::Object()),
        startCTime(0),
        startWTime(0),
        startBenchTime(0),
	ldbSteps(0)

{
    broadcast = new ControllerBroadcasts;
    reduction = ReductionMgr::Object()->willRequire(REDUCTIONS_BASIC);
    random = new Random(simParams->randomSeed);
    random->split(0,PatchMap::Object()->numPatches()+1);

    rescaleVelocities_sumTemps = 0;  rescaleVelocities_numTemps = 0;
    langevinPiston_strainRate = Tensor::diagonal(simParams->strainRate);
    if ( ! simParams->useFlexibleCell ) {
      BigReal avg = trace(langevinPiston_strainRate) / 3.;
      langevinPiston_strainRate = Tensor::identity(avg);
    }
}

Controller::~Controller(void)
{
    delete broadcast;
    delete reduction;
    delete random;
}

void Controller::threadRun(Controller* arg)
{
    arg->algorithm();
}

void Controller::run(void)
{
    // create a Thread and invoke it
    DebugM(4, "Starting thread in controller on this=" << this << "\n");
    thread = CthCreate((CthVoidFn)&(threadRun),(void*)(this),CTRL_STK_SZ);
    CthSetStrategyDefault(thread);
    awaken();
}

extern int eventEndOfTimeStep;

void Controller::algorithm(void)
{
  int scriptTask = SCRIPT_END;
  int scriptSeq = 0;
  if ( simParams->tclOn ) { Node::Object()->enableExitScheduler(); }
  while ( (! simParams->tclOn) ||
    (scriptTask = broadcast->scriptBarrier.get(scriptSeq++)) != SCRIPT_END) {
    switch ( scriptTask ) {
      case SCRIPT_RUN:
        break;
      case SCRIPT_MINIMIZE:
        minimize();
        continue;
      case SCRIPT_OUTPUT:
        enqueueCollections(FILE_OUTPUT);
        outputExtendedSystem(FILE_OUTPUT);
        Node::Object()->enableExitScheduler();
        continue;
      case SCRIPT_MEASURE:
        enqueueCollections(EVAL_MEASURE);
        Node::Object()->enableExitScheduler();
        continue;
      case SCRIPT_REINITVELS:
        iout << "REINITIALIZING VELOCITIES AT STEP " << simParams->firstTimestep
          << " TO " << simParams->initialTemp << " KELVIN.\n" << endi;
        Node::Object()->enableExitScheduler();
        continue;
      case SCRIPT_CHECKPOINT:
        iout << "CHECKPOINTING POSITIONS AT STEP " << simParams->firstTimestep
          << "\n" << endi;
        checkpoint_lattice = state->lattice;
        checkpoint_langevinPiston_strainRate = langevinPiston_strainRate;
        Node::Object()->enableExitScheduler();
        continue;
      case SCRIPT_REVERT:
        iout << "REVERTING POSITIONS AT STEP " << simParams->firstTimestep
          << "\n" << endi;
        state->lattice = checkpoint_lattice;
        langevinPiston_strainRate = checkpoint_langevinPiston_strainRate;
        Node::Object()->enableExitScheduler();
        continue;
    }

    if ( simParams->minimizeCGOn ) {
      minimize();
      if (! simParams->tclOn) break;
      Node::Object()->enableExitScheduler();
      continue;
    }

    int step = simParams->firstTimestep;

    const int numberOfSteps = simParams->N;
    const int stepsPerCycle = simParams->stepsPerCycle;

    nbondFreq = simParams->nonbondedFrequency;
    if ( simParams->fullDirectOn || simParams->FMAOn || simParams->PMEOn )
      slowFreq = simParams->fullElectFrequency;
    else
      slowFreq = simParams->nonbondedFrequency;

    receivePressure(step);
    reassignVelocities(step);
    printEnergies(step);
    traceUserEvent(eventEndOfTimeStep);
    outputExtendedSystem(step);
    rebalanceLoad(step);

    for ( ++step ; step <= numberOfSteps; ++step )
    {
        rescaleVelocities(step);
	tcoupleVelocities(step);
	berendsenPressure(step);
        enqueueCollections(step);
	langevinPiston1(step);
	receivePressure(step);
	langevinPiston2(step);
        reassignVelocities(step);
        printEnergies(step);
        traceUserEvent(eventEndOfTimeStep);
        outputExtendedSystem(step);
        cycleBarrier(!((step+1) % stepsPerCycle),step);
        rebalanceLoad(step);
    }

    if (! simParams->tclOn) break;
    Node::Object()->enableExitScheduler();
  }
  enqueueCollections(END_OF_RUN);
  outputExtendedSystem(END_OF_RUN);
  terminate();
}

#define CALCULATE \
  printMinimizeEnergies(step); \
  rebalanceLoad(step); \
  if ( step == numberOfSteps ) { \
    Node::Object()->enableExitScheduler(); \
    return; \
  } else ++step;

#define MOVETO(X) \
  if ( (X)-last.x ) { \
    if ( 0 ) { iout << "LINE MINIMIZER: MOVING FROM " << last.x << " TO " << (X) << "\n" << endi; } \
    broadcast->minimizeCoefficient.publish(minSeq++,(X)-last.x); \
    last.x = (X); \
    enqueueCollections(step); \
    CALCULATE \
    last.u = min_energy; \
    last.dudx = -1. * min_f_dot_v; \
  }

struct minpoint {
  BigReal x, u, dudx;
  minpoint() : x(0), u(0), dudx(0) { ; }
};

void Controller::minimize() {
  // iout << "Controller::minimize() called.\n" << endi;

  const int numberOfSteps = simParams->N;
  int step = simParams->firstTimestep;
  BigReal tinystep = simParams->minTinyStep;  // 1.0e-6
  BigReal babystep = simParams->minBabyStep;  // 1.0e-2
  BigReal linegoal = simParams->minLineGoal;  // 1.0e-4

  CALCULATE

  int minSeq = 0;
  int atStart = 2;
  BigReal old_f_dot_f = min_f_dot_f;
  broadcast->minimizeCoefficient.publish(minSeq++,0.);
  broadcast->minimizeCoefficient.publish(minSeq++,0.); // v = f
  min_f_dot_v = min_f_dot_f;
  min_v_dot_v = min_f_dot_f;
  int newDir = 1;
  while ( 1 ) {
    // line minimization
    // bracket minimum on line
    int bracket = 0;
    minpoint lo,hi,mid,last;
    BigReal x = last.x = 0;
    lo.x = x;
    lo.u = min_energy;
    lo.dudx = -1. * min_f_dot_v;
    mid = lo;
    BigReal tol = fabs( linegoal * min_f_dot_v );
    iout << "GRADIENT TOLERANCE: " << tol << "\n" << endi;
    x = babystep; if ( atStart ) { x = tinystep; }
    x *= sqrt( min_f_dot_f / min_v_dot_v ); MOVETO(x)
    // bracket minimum on line
    while ( last.u < mid.u ) {
      lo = mid; mid = last;
      x *= 2.0; MOVETO(x)
    }
    hi = last;
    iout << "BRACKET: " << (hi.x-lo.x) << " " << ((hi.u>lo.u?hi.u:lo.u)-mid.u) << " " << lo.dudx << " " << mid.dudx << " " << hi.dudx << " \n" << endi;
    // converge on minimum on line
    int itcnt = 0;
    while ( fabs(last.dudx) > tol && itcnt < 30 ) {
      // select new position
      if ( mid.dudx > 0. ) {
        BigReal altx = 0.1 * lo.x + 0.9 * mid.x;
        x = mid.dudx*(mid.x*mid.x-lo.x*lo.x) + 2*mid.x*(lo.u-mid.u);
        x /= 2*(mid.dudx*(mid.x-lo.x)+(lo.u-mid.u));
        if ( x > altx ) x = altx;
        if (x-last.x) { MOVETO(x) } else break;
        if ( last.u <= mid.u ) { hi = mid; mid = last; } else { lo = last; }
      } else {
        BigReal altx = 0.1 * hi.x + 0.9 * mid.x;
        x = mid.dudx*(mid.x*mid.x-hi.x*hi.x) + 2*mid.x*(hi.u-mid.u);
        x /= 2*(mid.dudx*(mid.x-hi.x)+(hi.u-mid.u));
        if ( x < altx ) x = altx;
        if (x-last.x) { MOVETO(x) } else break;
        if ( last.u <= mid.u ) { lo = mid; mid = last; } else { hi = last; }
      }
      iout << "BRACKET: " << (hi.x-lo.x) << " " << ((hi.u>lo.u?hi.u:lo.u)-mid.u) << " " << lo.dudx << " " << mid.dudx << " " << hi.dudx << " \n" << endi;
    }
    // new direction
    broadcast->minimizeCoefficient.publish(minSeq++,0.);
    BigReal c = min_f_dot_f / old_f_dot_f;
    iout << "NEW SEARCH DIRECTION " << c << "\n" << endi;
    c = ( c > 1.5 ? 1.5 : c );
    if ( atStart ) { c = 0; --atStart; }
    if ( c*c*min_v_dot_v > 100*min_f_dot_f ) {
      iout << "RESTARTING CONJUGATE GRADIENT ALGORITHM.\n" << endi;
      c = 0;
    }
    broadcast->minimizeCoefficient.publish(minSeq++,c); // v = c*v+f
    old_f_dot_f = min_f_dot_f;
    min_f_dot_v = c * min_f_dot_v + min_f_dot_f;
    min_v_dot_v = c*c*min_v_dot_v + 2*c*min_f_dot_v + min_f_dot_f;
  }

}

#undef MOVETO
#undef CALCULATE

void Controller::berendsenPressure(int step)
{
  const int freq = simParams->berendsenPressureFreq;
  if ( simParams->berendsenPressureOn && !((step-1)%freq) )
  {
    BigReal scalarPressure = trace(controlPressure) / 3.;
    BigReal factor = scalarPressure - simParams->berendsenPressureTarget;
    factor *= simParams->berendsenPressureCompressibility;
    factor *= ( simParams->dt * freq );
    factor /= simParams->berendsenPressureRelaxationTime;
    factor += 1.0;
    if ( factor < 0.9 ) {
      iout << iERROR << "Step " << step <<
	" volume rescaling factor limited to " <<
	0.9 << " from " << factor << "\n" << endi;
      factor = 0.9;
    }
    if ( factor > 1.1 ) {
      iout << iERROR << "Step " << step <<
	" volume rescaling factor limited to " <<
	1.1 << " from " << factor << "\n" << endi;
      factor = 1.1;
    }
    factor = cbrt(factor);
    broadcast->positionRescaleFactor.publish(step,Tensor::identity()*factor);
    state->lattice.rescale(Tensor::identity()*factor);
  }
}

void Controller::langevinPiston1(int step)
{
  if ( simParams->langevinPistonOn )
  {
    Tensor &strainRate = langevinPiston_strainRate;
    int cellDims = simParams->useFlexibleCell ? 1 : 3;
    BigReal dt = simParams->dt;
    BigReal dt_long = slowFreq * dt;
    BigReal kT = BOLTZMAN * simParams->langevinPistonTemp;
    BigReal tau = simParams->langevinPistonPeriod;
    BigReal mass = controlNumDegFreedom * kT * tau * tau * cellDims;

#ifdef DEBUG_PRESSURE
    iout << iINFO << "entering langevinPiston1, strain rate: " << strainRate << "\n";
#endif

    if ( ! ( (step-1) % slowFreq ) )
    {
      BigReal gamma = 1 / simParams->langevinPistonDecay;
      BigReal f1 = exp( -0.5 * dt_long * gamma );
      BigReal f2 = sqrt( ( 1. - f1*f1 ) * kT / mass );
      strainRate *= f1;
      if ( simParams->useFlexibleCell )
	strainRate += f2 * Tensor::diagonal(random->gaussian_vector());
      else
	strainRate += f2 * Tensor::identity(random->gaussian());
#ifdef DEBUG_PRESSURE
      iout << iINFO << "applying langevin, strain rate: " << strainRate << "\n";
#endif
    }

    // Apply surface tension.  If surfaceTensionTarget is zero, we get
    // the default (isotropic pressure) case.
    
    Tensor ptarget;
    ptarget.zz = simParams->langevinPistonTarget;
    ptarget.xx = ptarget.yy = simParams->langevinPistonTarget - 
        simParams->surfaceTensionTarget / state->lattice.c().z;

    strainRate += ( 0.5 * dt * cellDims * state->lattice.volume() / mass ) *
      ( controlPressure - ptarget );

#ifdef DEBUG_PRESSURE
    iout << iINFO << "integrating half step, strain rate: " << strainRate << "\n";
#endif

    if ( ! ( (step-1-slowFreq/2) % slowFreq ) )
    {
      // JCP FIX THIS:  NOT JUST ELEMENT-WISE EXP.  WHY NOT LINEAR?
      Tensor factor;
      factor.xx = exp( dt_long * strainRate.xx );
      factor.yy = exp( dt_long * strainRate.yy );
      factor.zz = exp( dt_long * strainRate.zz );
      broadcast->positionRescaleFactor.publish(step,factor);
      state->lattice.rescale(factor);
#ifdef DEBUG_PRESSURE
      iout << iINFO << "rescaling by: " << factor << "\n";
#endif
    }

    // corrections to integrator
    if ( ! ( step % nbondFreq ) )
    {
#ifdef DEBUG_PRESSURE
      iout << iINFO << "correcting strain rate for nbond, ";
#endif
      strainRate -= ( 0.5 * dt * cellDims * state->lattice.volume() / mass ) *
		( (nbondFreq - 1) * controlPressure_nbond );
#ifdef DEBUG_PRESSURE
      iout << "strain rate: " << strainRate << "\n";
#endif
    }
    if ( ! ( step % slowFreq ) )
    {
#ifdef DEBUG_PRESSURE
      iout << iINFO << "correcting strain rate for slow, ";
#endif
      strainRate -= ( 0.5 * dt * cellDims * state->lattice.volume() / mass ) *
		( (slowFreq - 1) * controlPressure_slow );
#ifdef DEBUG_PRESSURE
      iout << "strain rate: " << strainRate << "\n";
#endif
    }

  }
}

void Controller::langevinPiston2(int step)
{
  if ( simParams->langevinPistonOn )
  {
    Tensor &strainRate = langevinPiston_strainRate;
    int cellDims = simParams->useFlexibleCell ? 1 : 3;
    BigReal dt = simParams->dt;
    BigReal dt_long = slowFreq * dt;
    BigReal kT = BOLTZMAN * simParams->langevinPistonTemp;
    BigReal tau = simParams->langevinPistonPeriod;
    BigReal mass = controlNumDegFreedom * kT * tau * tau * cellDims;

    // corrections to integrator
    if ( ! ( step % nbondFreq ) )
    {
#ifdef DEBUG_PRESSURE
      iout << iINFO << "correcting strain rate for nbond, ";
#endif
      strainRate += ( 0.5 * dt * cellDims * state->lattice.volume() / mass ) *
		( (nbondFreq - 1) * controlPressure_nbond );
#ifdef DEBUG_PRESSURE
      iout << "strain rate: " << strainRate << "\n";
#endif
    }
    if ( ! ( step % slowFreq ) )
    {
#ifdef DEBUG_PRESSURE
      iout << iINFO << "correcting strain rate for slow, ";
#endif
      strainRate += ( 0.5 * dt * cellDims * state->lattice.volume() / mass ) *
		( (slowFreq - 1) * controlPressure_slow );
#ifdef DEBUG_PRESSURE
      iout << "strain rate: " << strainRate << "\n";
#endif
    }

    // Apply surface tension.  If surfaceTensionTarget is zero, we get
    // the default (isotropic pressure) case.
   
    Tensor ptarget;
    ptarget.zz = simParams->langevinPistonTarget;
    ptarget.xx = ptarget.yy = simParams->langevinPistonTarget -
        simParams->surfaceTensionTarget / state->lattice.c().z;

    strainRate += ( 0.5 * dt * cellDims * state->lattice.volume() / mass ) *
      ( controlPressure - ptarget );
 
#ifdef DEBUG_PRESSURE
    iout << iINFO << "integrating half step, strain rate: " << strainRate << "\n";
#endif

    if ( ! ( step % slowFreq ) )
    {
      BigReal gamma = 1 / simParams->langevinPistonDecay;
      BigReal f1 = exp( -0.5 * dt_long * gamma );
      BigReal f2 = sqrt( ( 1. - f1*f1 ) * kT / mass );
      strainRate *= f1;
      if ( simParams->useFlexibleCell )
	strainRate += f2 * Tensor::diagonal(random->gaussian_vector());
      else
	strainRate += f2 * Tensor::identity(random->gaussian());
#ifdef DEBUG_PRESSURE
      iout << iINFO << "applying langevin, strain rate: " << strainRate << "\n";
#endif
    }

#ifdef DEBUG_PRESSURE
    iout << iINFO << "exiting langevinPiston2, strain rate: " << strainRate << "\n";
#endif
  }
}

void Controller::rescaleVelocities(int step)
{
  const int rescaleFreq = simParams->rescaleFreq;
  if ( rescaleFreq > 0 ) {
    rescaleVelocities_sumTemps += temperature;  ++rescaleVelocities_numTemps;
    if ( rescaleVelocities_numTemps == rescaleFreq ) {
      BigReal avgTemp = rescaleVelocities_sumTemps / rescaleVelocities_numTemps;
      const BigReal rescaleTemp = simParams->rescaleTemp;
      BigReal factor = sqrt(rescaleTemp/avgTemp);
      broadcast->velocityRescaleFactor.publish(step,factor);
      iout << "RESCALING VELOCITIES AT STEP " << step
           << " FROM AVERAGE TEMPERATURE OF " << avgTemp
           << " TO " << rescaleTemp << " KELVIN.\n" << endi;
      rescaleVelocities_sumTemps = 0;  rescaleVelocities_numTemps = 0;
    }
  }
}

void Controller::reassignVelocities(int step)
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
    iout << "REASSIGNING VELOCITIES AT STEP " << step
         << " TO " << newTemp << " KELVIN.\n" << endi;
  }
}

void Controller::tcoupleVelocities(int step)
{
  if ( simParams->tCoupleOn )
  {
    const BigReal tCoupleTemp = simParams->tCoupleTemp;
    BigReal coefficient = 1.;
    if ( temperature > 0. ) coefficient = tCoupleTemp/temperature - 1.;
    broadcast->tcoupleCoefficient.publish(step,coefficient);
  }
}

static char *FORMAT(BigReal X)
{
  static char tmp_string[25];
  const double maxnum = 1000000000000.0;
  if ( X > maxnum ) X = maxnum;
  if ( X < -maxnum ) X = -maxnum;
  sprintf(tmp_string,"%.4f ",X); 
  NAMD_pad(tmp_string, 12);
  return  tmp_string;
}

static char *ETITLE(int X)
{
  static char tmp_string[21];
  sprintf(tmp_string,"ENERGY: %6d ",X); 
  return  tmp_string;
}

void Controller::receivePressure(int step)
{
    Node *node = Node::Object();
    Molecule *molecule = node->molecule;
    SimParameters *simParameters = node->simParameters;
    Lattice &lattice = state->lattice;

    reduction->require();

    // BigReal intKineticEnergy;
    Tensor virial;
    Tensor virial_normal;
    Tensor virial_nbond;
    Tensor virial_slow;
#ifdef ALTVIRIAL
    Tensor altVirial_normal;
    Tensor altVirial_nbond;
    Tensor altVirial_slow;
#endif
    Tensor intVirial;
    Tensor intVirial_normal;
    Tensor intVirial_nbond;
    Tensor intVirial_slow;
    BigReal volume;

    int numAtoms = molecule->numAtoms;
    numDegFreedom = 3 * numAtoms;
    int numFixedAtoms = molecule->numFixedAtoms;
    if ( numFixedAtoms ) numDegFreedom -= 3 * numFixedAtoms;
    if ( ! ( numFixedAtoms || molecule->numConstraints
	|| simParameters->comMove || simParameters->langevinOn ) ) {
      numDegFreedom -= 3;
    }
    int numRigidBonds = molecule->numRigidBonds;
    int numFixedRigidBonds = molecule->numFixedRigidBonds;
    numDegFreedom -= ( numRigidBonds - numFixedRigidBonds );

    kineticEnergy = reduction->item(REDUCTION_KINETIC_ENERGY);
    // intKineticEnergy = reduction->item(REDUCTION_INT_KINETIC_ENERGY);

    GET_TENSOR(virial_normal,reduction,REDUCTION_VIRIAL_NORMAL);
    GET_TENSOR(virial_nbond,reduction,REDUCTION_VIRIAL_NBOND);
    GET_TENSOR(virial_slow,reduction,REDUCTION_VIRIAL_SLOW);

#ifdef ALTVIRIAL
    GET_TENSOR(altVirial_normal,reduction,REDUCTION_ALT_VIRIAL_NORMAL);
    GET_TENSOR(altVirial_nbond,reduction,REDUCTION_ALT_VIRIAL_NBOND);
    GET_TENSOR(altVirial_slow,reduction,REDUCTION_ALT_VIRIAL_SLOW);
#endif

    GET_TENSOR(intVirial_normal,reduction,REDUCTION_INT_VIRIAL_NORMAL);
    GET_TENSOR(intVirial_nbond,reduction,REDUCTION_INT_VIRIAL_NBOND);
    GET_TENSOR(intVirial_slow,reduction,REDUCTION_INT_VIRIAL_SLOW);

    temperature = 2.0 * kineticEnergy / ( numDegFreedom * BOLTZMAN );

    if ( (volume=lattice.volume()) != 0. )
    {
      // kinetic energy component included in virials
      pressure_normal = virial_normal / volume;
      groupPressure_normal = ( virial_normal - intVirial_normal ) / volume;

      if ( ! ( step % nbondFreq ) )
      {
        pressure_nbond = virial_nbond / volume;
        groupPressure_nbond = ( virial_nbond - intVirial_nbond ) / volume;
      }

      if ( ! ( step % slowFreq ) )
      {
        pressure_slow = virial_slow / volume;
        groupPressure_slow = ( virial_slow - intVirial_slow ) / volume;
      }

/*
      iout << "VIRIALS: " << virial_normal << " " << virial_nbond << " " <<
	virial_slow << " " << ( virial_normal - intVirial_normal ) << " " <<
	( virial_nbond - intVirial_nbond ) << " " <<
	( virial_slow - intVirial_slow ) << "\n";
*/

      pressure = pressure_normal + pressure_nbond + pressure_slow; 
      groupPressure = groupPressure_normal + groupPressure_nbond +
						groupPressure_slow;
    }
    else
    {
      pressure = Tensor();
      groupPressure = Tensor();
    }

    if ( simParameters->useGroupPressure )
    {
      controlPressure_normal = groupPressure_normal;
      controlPressure_nbond = groupPressure_nbond;
      controlPressure_slow = groupPressure_slow;
      controlPressure = groupPressure;
      controlNumDegFreedom = molecule->numHydrogenGroups;
      if ( ! ( numFixedAtoms || molecule->numConstraints
	|| simParameters->comMove || simParameters->langevinOn ) ) {
        controlNumDegFreedom -= 1;
      }
    }
    else
    {
      controlPressure_normal = pressure_normal;
      controlPressure_nbond = pressure_nbond;
      controlPressure_slow = pressure_slow;
      controlPressure = pressure;
      controlNumDegFreedom = numDegFreedom / 3;
    }

    if ( ! simParameters->useFlexibleCell ) {
      controlPressure_normal =
		Tensor::identity(trace(controlPressure_normal)/3.);
      controlPressure_nbond =
		Tensor::identity(trace(controlPressure_nbond)/3.);
      controlPressure_slow =
		Tensor::identity(trace(controlPressure_slow)/3.);
      controlPressure =
		Tensor::identity(trace(controlPressure)/3.);
    }

#ifdef DEBUG_PRESSURE
    iout << iINFO << "Control pressure = " << controlPressure <<
      " = " << controlPressure_normal << " + " <<
      controlPressure_nbond << " + " << controlPressure_slow << "\n" << endi;
#endif

}

void Controller::compareChecksums(int step) {
    Node *node = Node::Object();
    Molecule *molecule = node->molecule;

    // Some basic correctness checking
    BigReal checksum;

    checksum = reduction->item(REDUCTION_ATOM_CHECKSUM);
    if ( ((int)checksum) != molecule->numAtoms )
      NAMD_bug("Bad global atom count!\n");

    checksum = reduction->item(REDUCTION_COMPUTE_CHECKSUM);
    if ( ((int)checksum) != computeChecksum ) {
      if ( computeChecksum )
        NAMD_bug("Bad global compute count!\n");
      else
        computeChecksum = ((int)checksum);
    }

    checksum = reduction->item(REDUCTION_BOND_CHECKSUM);
    if ( ((int)checksum) != molecule->numCalcBonds )
      NAMD_bug("Bad global bond count!\n");

    checksum = reduction->item(REDUCTION_ANGLE_CHECKSUM);
    if ( ((int)checksum) != molecule->numCalcAngles )
      NAMD_bug("Bad global angle count!\n");

    checksum = reduction->item(REDUCTION_DIHEDRAL_CHECKSUM);
    if ( ((int)checksum) != molecule->numCalcDihedrals )
      NAMD_bug("Bad global dihedral count!\n");

    checksum = reduction->item(REDUCTION_IMPROPER_CHECKSUM);
    if ( ((int)checksum) != molecule->numCalcImpropers )
      NAMD_bug("Bad global improper count!\n");

/*
    checksum = reduction->item(REDUCTION_EXCLUSION_CHECKSUM);
    if ( ((int)checksum) != molecule->numCalcExclusions )
      NAMD_bug("Bad global exclusion count!\n");
*/

    checksum = reduction->item(REDUCTION_MARGIN_VIOLATIONS);
    if ( ((int)checksum) ) iout << iWARN << ((int)checksum) <<
        " margin violations detected during timestep " << step << ".\n" << endi;
}

void Controller::printMinimizeEnergies(int step) {
    reduction->require();
    compareChecksums(step);

    BigReal bondEnergy;
    BigReal angleEnergy;
    BigReal dihedralEnergy;
    BigReal improperEnergy;
    BigReal boundaryEnergy;
    BigReal miscEnergy;
    BigReal totalEnergy;

    bondEnergy = reduction->item(REDUCTION_BOND_ENERGY);
    angleEnergy = reduction->item(REDUCTION_ANGLE_ENERGY);
    dihedralEnergy = reduction->item(REDUCTION_DIHEDRAL_ENERGY);
    improperEnergy = reduction->item(REDUCTION_IMPROPER_ENERGY);
    boundaryEnergy = reduction->item(REDUCTION_BC_ENERGY);
    miscEnergy = reduction->item(REDUCTION_MISC_ENERGY);
    electEnergy = reduction->item(REDUCTION_ELECT_ENERGY);
    ljEnergy = reduction->item(REDUCTION_LJ_ENERGY);
    electEnergySlow = reduction->item(REDUCTION_ELECT_ENERGY_SLOW);

    totalEnergy = bondEnergy + angleEnergy + dihedralEnergy + improperEnergy +
	electEnergy + electEnergySlow + ljEnergy + boundaryEnergy + miscEnergy;

    if ( ( step % 10 ) == 0 ) {
	iout << "ETITLE:     TS    BOND        ANGLE       "
	     << "DIHED       IMPRP       ELECT       VDW       "
	     << "BOUNDARY    MISC        TOTAL\n" << endi;
    }

    // N.B.  HP's aCC compiler merges FORMAT calls in the same expression.
    //       Need separate statements because data returned in static array.

    iout << ETITLE(step);
    iout << FORMAT(bondEnergy);
    iout << FORMAT(angleEnergy);
    iout << FORMAT(dihedralEnergy);
    iout << FORMAT(improperEnergy);
    iout << FORMAT(electEnergy+electEnergySlow);
    iout << FORMAT(ljEnergy);
    iout << FORMAT(boundaryEnergy);
    iout << FORMAT(miscEnergy);
    iout << FORMAT(totalEnergy);
    iout << "\n" << endi;

    min_energy = totalEnergy;
    min_f_dot_f = reduction->item(REDUCTION_MIN_F_DOT_F);
    min_f_dot_v = reduction->item(REDUCTION_MIN_F_DOT_V);
    min_v_dot_v = reduction->item(REDUCTION_MIN_V_DOT_V);
}

void Controller::printEnergies(int step)
{
    Node *node = Node::Object();
    Molecule *molecule = node->molecule;
    SimParameters *simParameters = node->simParameters;
    Lattice &lattice = state->lattice;

    compareChecksums(step);

    BigReal bondEnergy;
    BigReal angleEnergy;
    BigReal dihedralEnergy;
    BigReal improperEnergy;
    BigReal boundaryEnergy;
    BigReal miscEnergy;
    BigReal totalEnergy;
    Vector momentum;
    Vector angularMomentum;
    BigReal volume = lattice.volume();

    bondEnergy = reduction->item(REDUCTION_BOND_ENERGY);
    angleEnergy = reduction->item(REDUCTION_ANGLE_ENERGY);
    dihedralEnergy = reduction->item(REDUCTION_DIHEDRAL_ENERGY);
    improperEnergy = reduction->item(REDUCTION_IMPROPER_ENERGY);
    boundaryEnergy = reduction->item(REDUCTION_BC_ENERGY);
    miscEnergy = reduction->item(REDUCTION_MISC_ENERGY);

    if ( ! ( step % nbondFreq ) )
    {
      electEnergy = reduction->item(REDUCTION_ELECT_ENERGY);
      ljEnergy = reduction->item(REDUCTION_LJ_ENERGY);
    }

    if ( ! ( step % slowFreq ) )
    {
      electEnergySlow = reduction->item(REDUCTION_ELECT_ENERGY_SLOW);
    }

    momentum.x = reduction->item(REDUCTION_MOMENTUM_X);
    momentum.y = reduction->item(REDUCTION_MOMENTUM_Y);
    momentum.z = reduction->item(REDUCTION_MOMENTUM_Z);
    angularMomentum.x = reduction->item(REDUCTION_ANGULAR_MOMENTUM_X);
    angularMomentum.y = reduction->item(REDUCTION_ANGULAR_MOMENTUM_Y);
    angularMomentum.z = reduction->item(REDUCTION_ANGULAR_MOMENTUM_Z);

    totalEnergy = bondEnergy + angleEnergy + dihedralEnergy + improperEnergy +
	electEnergy + electEnergySlow + ljEnergy + kineticEnergy +
	boundaryEnergy + miscEnergy;

    if ( node->simParameters->outputMomenta &&
         ! ( step % node->simParameters->outputMomenta ) )
    {
      iout << "MOMENTUM: " << step 
           << " P: " << momentum
           << " L: " << angularMomentum
           << "\n" << endi;
    }

    if ( node->simParameters->outputPressure &&
         ! ( step % node->simParameters->outputPressure ) )
    {
      iout << "PRESSURE: " << step << " "
           << PRESSUREFACTOR * pressure << "\n"
           << "GPRESSURE: " << step << " "
           << PRESSUREFACTOR * groupPressure << "\n" << endi;
    }

    if (node->simParameters->IMDon && !(step % node->simParameters->IMDfreq)) {
      IMDEnergies energies;
      energies.tstep = step;
      energies.T = temperature;
      energies.Etot = totalEnergy;
      energies.Epot = totalEnergy - kineticEnergy;
      energies.Evdw = ljEnergy;
      energies.Eelec = electEnergy + electEnergySlow;
      energies.Ebond = bondEnergy;
      energies.Eangle = angleEnergy;
      energies.Edihe = dihedralEnergy;
      energies.Eimpr = improperEnergy;
      Node::Object()->imd->gather_energies(&energies);
    }
  
    int stepInRun = step - simParams->firstTimestep;
    int benchPhase;
    if ( stepInRun % simParams->firstLdbStep == 0 )
    switch ( benchPhase = stepInRun / simParams->firstLdbStep )
    {
    case 0:
    case 2:
    case 4:
      startBenchTime = CmiWallTimer();
      break;
    case 1:
    case 3:
    case 5:
      iout << iINFO;
      if ( benchPhase == 1 ) iout << "Initial time: ";
      else iout << "Benchmark time: ";
      {
        BigReal wallTime = CmiWallTimer() - startBenchTime;
        BigReal wallPerStep =
		(CmiWallTimer() - startBenchTime) / simParams->firstLdbStep;
	BigReal ns = simParams->dt / 1000000.0;
	BigReal days = 1.0 / (24.0 * 60.0 * 60.0);
	BigReal daysPerNano = wallPerStep * days / ns;
	iout << wallPerStep << " s/step " << daysPerNano << " days/ns ";
        iout << (memusage()/1024) << " kB of memory in use.\n" << endi;
      }
      break;
    }

    if ( simParams->outputTiming && ! ( step % simParams->outputTiming ) )
    {
      const double endWTime = CmiWallTimer();
      const double endCTime = CmiTimer();

      const double elapsedW = 
	(endWTime - startWTime) / simParams->outputTiming;
      const double elapsedC = 
	(endCTime - startCTime) / simParams->outputTiming;

      const double remainingW = elapsedW * (simParams->N - step);
      const double remainingW_hours = remainingW / 3600;

      startWTime = endWTime;
      startCTime = endCTime;

      if ( step >= (simParams->firstTimestep + simParams->outputTiming) ) {
        iout << "TIMING: " << step
             << "  CPU: " << endCTime << ", " << elapsedC << "/step"
             << "  Wall: " << endWTime << ", " << elapsedW << "/step"
             << ", " << remainingW_hours << " hours remaining"
             << ", " << (memusage()/1024) << " kB of memory in use"
             << ".\n" << endi;
      }
    }

    // callback to Tcl with whatever we can
#ifdef NAMD_TCL
#define CALLBACKDATA(LABEL,VALUE) \
		labels << (LABEL) << " "; values << (VALUE) << " ";
#define CALLBACKLIST(LABEL,VALUE) \
		labels << (LABEL) << " "; values << "{" << (VALUE) << "} ";
    if (node->getScript() && node->getScript()->doCallback()) {
      ostrstream labels, values;
      CALLBACKDATA("TS",step);
      CALLBACKDATA("BOND",bondEnergy);
      CALLBACKDATA("ANGLE",angleEnergy);
      CALLBACKDATA("DIHED",dihedralEnergy);
      CALLBACKDATA("IMPRP",improperEnergy);
      CALLBACKDATA("ELECT",electEnergy+electEnergySlow);
      CALLBACKDATA("VDW",ljEnergy);
      CALLBACKDATA("BOUNDARY",boundaryEnergy);
      CALLBACKDATA("MISC",miscEnergy);
      CALLBACKDATA("KINETIC",kineticEnergy);
      CALLBACKDATA("TOTAL",totalEnergy);
      CALLBACKDATA("TEMP",temperature);
      CALLBACKLIST("PRESSURE",pressure*PRESSUREFACTOR);
      CALLBACKLIST("GPRESSURE",groupPressure*PRESSUREFACTOR);
      CALLBACKDATA("VOLUME",lattice.volume());
      CALLBACKLIST("CELL_A",lattice.a());
      CALLBACKLIST("CELL_B",lattice.b());
      CALLBACKLIST("CELL_C",lattice.c());
      CALLBACKLIST("CELL_O",lattice.origin());
      labels << "PERIODIC"; values << "{" << lattice.a_p() << " "
		<< lattice.b_p() << " " << lattice.c_p() << "}";

      char *labelstr = labels.str();
      char *valuestr = values.str();
      node->getScript()->doCallback(labelstr,valuestr);
      delete [] labelstr;
      delete [] valuestr;
    }
#undef CALLBACKDATA
#endif

    // NO CALCULATIONS OR REDUCTIONS BEYOND THIS POINT!!!
    if ( step % node->simParameters->outputEnergies ) return;
    // ONLY OUTPUT SHOULD OCCUR BELOW THIS LINE!!!

    int printAtomicPressure = 1;
#ifndef DEBUG_PRESSURE
    // if ( simParams->rigidBonds != RIGID_NONE ) { printAtomicPressure = 0; }
#endif

    if ( (step % (10 * node->simParameters->outputEnergies) ) == 0 )
    {
	iout << "ETITLE:     TS    BOND        ANGLE       "
	     << "DIHED       IMPRP       ELECT       VDW       "
	     << "BOUNDARY    MISC        KINETIC        TOTAL     TEMP";
	if ( volume != 0. ) {
	  if ( printAtomicPressure ) iout << "     PRESSURE";
	  iout << "    GPRESSURE";
	  iout << "    VOLUME";
	}
	iout << "\n" << endi;
    }

    // N.B.  HP's aCC compiler merges FORMAT calls in the same expression.
    //       Need separate statements because data returned in static array.

    iout << ETITLE(step);
    iout << FORMAT(bondEnergy);
    iout << FORMAT(angleEnergy);
    iout << FORMAT(dihedralEnergy);
    iout << FORMAT(improperEnergy);
    iout << FORMAT(electEnergy+electEnergySlow);
    iout << FORMAT(ljEnergy);
    iout << FORMAT(boundaryEnergy);
    iout << FORMAT(miscEnergy);
    iout << FORMAT(kineticEnergy);
    iout << FORMAT(totalEnergy);
    iout << FORMAT(temperature);

#if(CMK_CCS_AVAILABLE)
     char webout[80];
     sprintf(webout,"%d %d %d %d",(int)totalEnergy,
	     (int)(totalEnergy - kineticEnergy),
	     (int)kineticEnergy,(int)temperature);
     CApplicationDepositNode0Data(webout);
#endif

    if ( volume != 0. )
    {
	if ( printAtomicPressure ) {
	  iout << FORMAT(trace(pressure)*PRESSUREFACTOR/3.);
	}
	iout << FORMAT(trace(groupPressure)*PRESSUREFACTOR/3.);
	iout << FORMAT(volume);
    }

    iout << "\n" << endi;
}

void Controller::writeExtendedSystemLabels(ofstream &file) {
  Lattice &lattice = state->lattice;
  file << "#$LABELS step";
  if ( lattice.a_p() ) file << " a_x a_y a_z";
  if ( lattice.b_p() ) file << " b_x b_y b_z";
  if ( lattice.c_p() ) file << " c_x c_y c_z";
  file << " o_x o_y o_z";
  if ( simParams->langevinPistonOn ) {
    file << " s_x s_y s_z s_u s_v s_w";
  }
  file << endl;
}

void Controller::writeExtendedSystemData(int step, ofstream &file) {
  Lattice &lattice = state->lattice;
  file << step;
    if ( lattice.a_p() ) file << " " << lattice.a().x << " " << lattice.a().y << " " << lattice.a().z;
    if ( lattice.b_p() ) file << " " << lattice.b().x << " " << lattice.b().y << " " << lattice.b().z;
    if ( lattice.c_p() ) file << " " << lattice.c().x << " " << lattice.c().y << " " << lattice.c().z;
    file << " " << lattice.origin().x << " " << lattice.origin().y << " " << lattice.origin().z;
  if ( simParams->langevinPistonOn ) {
    Vector strainRate = diagonal(langevinPiston_strainRate);
    Vector strainRate2 = off_diagonal(langevinPiston_strainRate);
    file << " " << strainRate.x;
    file << " " << strainRate.y;
    file << " " << strainRate.z;
    file << " " << strainRate2.x;
    file << " " << strainRate2.y;
    file << " " << strainRate2.z;
  }
  file << endl;
}

void Controller::enqueueCollections(int timestep)
{
  if ( Output::coordinateNeeded(timestep) )
    collection->enqueuePositions(timestep);
  if ( Output::velocityNeeded(timestep) )
    collection->enqueueVelocities(timestep);
}

void Controller::outputExtendedSystem(int step)
{

  if ( step >= 0 ) {

    // Write out eXtended System Trajectory (XST) file
    if ( simParams->xstFrequency &&
         ((step % simParams->xstFrequency) == 0) )
    {
      if ( ! xstFile.rdbuf()->is_open() )
      {
        NAMD_backup_file(simParams->xstFilename);
        xstFile.open(simParams->xstFilename);
        iout << "OPENING EXTENDED SYSTEM TRAJECTORY FILE\n" << endi;
        xstFile << "# NAMD extended system trajectory file" << endl;
        writeExtendedSystemLabels(xstFile);
      }
      writeExtendedSystemData(step,xstFile);
      xstFile.flush();
    }

    // Write out eXtended System Configuration (XSC) files
    //  Output a restart file
    if ( simParams->restartFrequency &&
         ((step % simParams->restartFrequency) == 0) &&
         (step != simParams->firstTimestep) )
    {
      char fname[140];
      strcpy(fname, simParams->restartFilename);
      strcat(fname, ".xsc");
      NAMD_backup_file(fname,".old");
      ofstream xscFile(fname);
      iout << "WRITING EXTENDED SYSTEM TO RESTART FILE AT STEP "
		<< step << "\n" << endi;
      xscFile << "# NAMD extended system configuration restart file" << endl;
      writeExtendedSystemLabels(xscFile);
      writeExtendedSystemData(step,xscFile);
    }

  }

  //  Output final coordinates
  if (step == FILE_OUTPUT || step == END_OF_RUN)
  {
    static char fname[140];
    strcpy(fname, simParams->outputFilename);
    strcat(fname, ".xsc");
    NAMD_backup_file(fname);
    ofstream xscFile(fname);
    iout << "WRITING EXTENDED SYSTEM TO OUTPUT FILE AT STEP "
		<< simParams->N << "\n" << endi;
    xscFile << "# NAMD extended system configuration output file" << endl;
    writeExtendedSystemLabels(xscFile);
    writeExtendedSystemData(simParams->N,xscFile);
  }

  //  Close trajectory file
  if (step == END_OF_RUN) {
    if ( xstFile.rdbuf()->is_open() ) {
      xstFile.close();
      iout << "CLOSING EXTENDED SYSTEM TRAJECTORY FILE\n" << endi;
    }
  }

}

void Controller::rebalanceLoad(int)
{
  if ( ! ldbSteps ) { 
    ldbSteps = LdbCoordinator::Object()->steps();
  }
  if ( ! --ldbSteps ) {
    LdbCoordinator::Object()->rebalance(this);
  }
}

void Controller::cycleBarrier(int doBarrier, int step) {
#ifdef CYCLE_BARRIER
	if (doBarrier) {
	  broadcast->cycleBarrier.publish(step,1);
	  CkPrintf("Cycle time at sync Wall: %f CPU %f\n",
		  CmiWallTimer(),CmiTimer());
	}
#endif
}

void Controller::terminate(void) {
  Node::Object()->enableHaltBarrier();
  CthFree(thread);
  CthSuspend();
}

