/**
***  Copyright (c) 1995, 1996, 1997, 1998, 1999, 2000 by
***  The Board of Trustees of the University of Illinois.
***  All rights reserved.
**/

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

#ifdef NAMD_CCS
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
        startBenchTime(0)

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
    arg->algorithm(0);
}

void Controller::run(void)
{
    // create a Thread and invoke it
    DebugM(4, "Starting thread in controller on this=" << this << "\n");
#ifdef NAMD_TCL
    thread = Node::Object()->getScript()->thread;
#else
    thread = CthCreate((CthVoidFn)&(threadRun),(void*)(this),CTRL_STK_SZ);
    CthSetStrategyDefault(thread);
#endif
    awaken();
}

extern int eventEndOfTimeStep;

void Controller::algorithm(int task)
{
  if (simParams->tclOn) broadcast->scriptBarrier.publish(scriptSeq++,task);

  switch ( task ) {
    case 0:
      if ( simParams->tclOn ) {
        enqueueCollections(0);
        outputExtendedSystem(-1);
        return;
      }
    case 1:
      break;
    case 2:
      collection->enqueuePositions(0);
      collection->enqueueVelocities(0);
      outputExtendedSystem(-1);
      return;
    }

    int step = simParams->firstTimestep;
    int first = 1;

    const int numberOfSteps = simParams->N;

    nbondFreq = simParams->nonbondedFrequency;
    if ( simParams->fullDirectOn || simParams->FMAOn || simParams->PMEOn )
      slowFreq = simParams->fullElectFrequency;
    else
      slowFreq = simParams->nonbondedFrequency;

    for ( ; step <= numberOfSteps; ++step, first = 0 )
    {
        if ( ! first ) rescaleVelocities(step);
	if ( ! first ) tcoupleVelocities(step);
	if ( ! first ) berendsenPressure(step);
        if ( ! first ) enqueueCollections(step);
        traceUserEvent(eventEndOfTimeStep);
	if ( ! first ) langevinPiston1(step);
	receivePressure(step);
	if ( ! first ) langevinPiston2(step);
        reassignVelocities(step);
        printEnergies(step);
        outputExtendedSystem(step);
        //rescaleVelocities(step);
	//tcoupleVelocities(step);
	//berendsenPressure(step);
#ifdef CYCLE_BARRIER
	if (!((step+1) % stepsPerCycle))
	{
	  broadcast->cycleBarrier.publish(step,1);
	  CkPrintf("Cycle time at sync Wall: %f CPU %f\n",
		  CmiWallTimer(),CmiTimer());
	}
#endif
	if ( LdbCoordinator::Object()->balanceNow(step) ) {
	  LdbCoordinator::Object()->rebalance(this);
	}
    }

  if ( ! task ) {
    enqueueCollections(0);
    outputExtendedSystem(-1);
    terminate();
  }
}

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

    strainRate += ( 0.5 * dt * cellDims * state->lattice.volume() / mass ) *
	( controlPressure - Tensor::identity(simParams->langevinPistonTarget) );

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

    strainRate += ( 0.5 * dt * cellDims * state->lattice.volume() / mass ) *
	( controlPressure - Tensor::identity(simParams->langevinPistonTarget) );

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
  static char tmp_string[21];
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
    Tensor altVirial_normal;
    Tensor altVirial_nbond;
    Tensor altVirial_slow;
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

    GET_TENSOR(altVirial_normal,reduction,REDUCTION_ALT_VIRIAL_NORMAL);
    GET_TENSOR(altVirial_nbond,reduction,REDUCTION_ALT_VIRIAL_NBOND);
    GET_TENSOR(altVirial_slow,reduction,REDUCTION_ALT_VIRIAL_SLOW);

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

void Controller::printEnergies(int step)
{
    Node *node = Node::Object();
    Molecule *molecule = node->molecule;
    SimParameters *simParameters = node->simParameters;
    Lattice &lattice = state->lattice;

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

    checksum = reduction->item(REDUCTION_EXCLUSION_CHECKSUM);
    if ( ((int)checksum) != molecule->numCalcExclusions )
      NAMD_bug("Bad global exclusion count!\n");

    checksum = reduction->item(REDUCTION_MARGIN_VIOLATIONS);
    if ( ((int)checksum) ) iout << iWARN << ((int)checksum) <<
        " margin violations detected during timestep " << step << ".\n" << endi;

    BigReal bondEnergy;
    BigReal angleEnergy;
    BigReal dihedralEnergy;
    BigReal improperEnergy;
    BigReal electEnergy;
    BigReal electEnergySlow;
    BigReal ljEnergy;
    BigReal boundaryEnergy;
    BigReal miscEnergy;
    BigReal smdEnergy;
    BigReal totalEnergy;
    Vector momentum;
    Vector angularMomentum;
    BigReal volume = lattice.volume();

    bondEnergy = reduction->item(REDUCTION_BOND_ENERGY);
    angleEnergy = reduction->item(REDUCTION_ANGLE_ENERGY);
    dihedralEnergy = reduction->item(REDUCTION_DIHEDRAL_ENERGY);
    improperEnergy = reduction->item(REDUCTION_IMPROPER_ENERGY);
    electEnergy = reduction->item(REDUCTION_ELECT_ENERGY);
    electEnergySlow = reduction->item(REDUCTION_ELECT_ENERGY_SLOW);
    ljEnergy = reduction->item(REDUCTION_LJ_ENERGY);
    boundaryEnergy = reduction->item(REDUCTION_BC_ENERGY);
    miscEnergy = reduction->item(REDUCTION_MISC_ENERGY);
    smdEnergy = reduction->item(REDUCTION_SMD_ENERGY);

    momentum.x = reduction->item(REDUCTION_MOMENTUM_X);
    momentum.y = reduction->item(REDUCTION_MOMENTUM_Y);
    momentum.z = reduction->item(REDUCTION_MOMENTUM_Z);
    angularMomentum.x = reduction->item(REDUCTION_ANGULAR_MOMENTUM_X);
    angularMomentum.y = reduction->item(REDUCTION_ANGULAR_MOMENTUM_Y);
    angularMomentum.z = reduction->item(REDUCTION_ANGULAR_MOMENTUM_Z);

    totalEnergy = bondEnergy + angleEnergy + dihedralEnergy + improperEnergy +
	electEnergy + electEnergySlow + ljEnergy + kineticEnergy +
	boundaryEnergy + miscEnergy + smdEnergy;

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
      energies.T = temperature;
      energies.Etot = totalEnergy;
      energies.Epot = totalEnergy - kineticEnergy;
      energies.Evdw = ljEnergy;
      energies.Eelec = electEnergy + electEnergySlow;
      energies.Ebond = bondEnergy;
      energies.Eangle = angleEnergy;
      energies.Edihe = dihedralEnergy;
      energies.Eimpr = improperEnergy;
      Node::Object()->imd->gather_energies(step, &energies);
    }
  
    int stepInRun = step - simParams->firstTimestep;
    int benchPhase;
    if ( stepInRun % simParams->firstLdbStep == 0 )
    switch ( benchPhase = stepInRun / simParams->firstLdbStep )
    {
    case 0:
    case 2:
      startBenchTime = CmiWallTimer();
      break;
    case 1:
    case 3:
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
	iout << wallPerStep << " s/step "
		<< daysPerNano << " days/ns\n" << endi;
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
             << ", " << remainingW_hours << " hours remaining\n" << endi;
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

      const char *labelstr = labels.str();
      const char *valuestr = values.str();
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
	if (simParams->SMDOn) iout << "     SMD";
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

#ifdef NAMD_CCS
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

    if (simParams->SMDOn) {
      iout << FORMAT(smdEnergy);
    }

    iout << "\n" << endi;
}

void Controller::writeExtendedSystemLabels(ofstream &file) {
  file << "#$LABELS step a_x a_y a_z b_x b_y b_z c_x c_y c_z o_x o_y o_z";
  if ( simParams->langevinPistonOn ) {
    file << " s_x s_y s_z s_u s_v s_w";
  }
  file << endl;
}

void Controller::writeExtendedSystemData(int step, ofstream &file) {
  Lattice &lattice = state->lattice;
  file << step
    << " " << lattice.a().x << " " << lattice.a().y << " " << lattice.a().z
    << " " << lattice.b().x << " " << lattice.b().y << " " << lattice.b().z
    << " " << lattice.c().x << " " << lattice.c().y << " " << lattice.c().z
    << " " << lattice.origin().x << " " << lattice.origin().y << " " << lattice.origin().z;
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
    // Write out eXtended System Trajectory (XST) file
    if ( step != -1 && simParams->xstFrequency != -1 &&
         ! ( step % simParams->xstFrequency ) )
    {
      if ( step == simParams->firstTimestep )
      {
        NAMD_backup_file(simParams->xstFilename);
        xstFile.open(simParams->xstFilename);
        xstFile << "# NAMD extended system trajectory file" << endl;
        writeExtendedSystemLabels(xstFile);
      }
      writeExtendedSystemData(step,xstFile);
      xstFile.flush();
      if ( simParams->xstFrequency != -1 && step == simParams->N )
      {
        xstFile.close();
      }
    }

    // Write out eXtended System Configuration (XSC) files
    //  Output a restart file
    if ( step != -1 && (simParams->restartFrequency != -1) &&
         ((step % simParams->restartFrequency) == 0) &&
         (step != simParams->firstTimestep) )
    {
      char fname[140];
      strcpy(fname, simParams->restartFilename);
      strcat(fname, ".xsc");
      NAMD_backup_file(fname);
      ofstream xscFile(fname);
      iout << "WRITING EXTENDED SYSTEM TO RESTART FILE AT STEP "
		<< step << "\n" << endi;
      xscFile << "# NAMD extended system configuration restart file" << endl;
      writeExtendedSystemLabels(xscFile);
      writeExtendedSystemData(step,xscFile);
    }

    //  Output final coordinates
    if (step == -1)
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
}

void Controller::terminate(void) {
  Node::Object()->enableHaltBarrier();
  CthFree(thread);
  CthSuspend();
}

