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

#ifdef NAMDCCS
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
    langevinPiston_strainRate = simParams->strainRate;
    if ( ! simParams->useFlexibleCell ) {
      BigReal avg = (langevinPiston_strainRate * Vector(1,1,1)) / 3.;
      langevinPiston_strainRate = Vector(avg,avg,avg);
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
  int scriptTask = 1;
  int scriptSeq = 0;
  if (simParams->tclOn) Node::Object()->enableScriptBarrier();
  while ( (! simParams->tclOn) ||
	(scriptTask = broadcast->scriptBarrier.get(scriptSeq++)) ) {
    switch ( scriptTask ) {
      case 1:
        break;
      case 2:
	collection->enqueuePositions(0);
	Node::Object()->enableScriptBarrier();
	continue;
      case 3:
	collection->enqueueVelocities(0);
	Node::Object()->enableScriptBarrier();
	continue;
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

    if (simParams->tclOn) Node::Object()->enableScriptBarrier();
    else break;
  }

  enqueueCollections(0);
  terminate();
}

void Controller::berendsenPressure(int step)
{
  const int freq = simParams->berendsenPressureFreq;
  if ( simParams->berendsenPressureOn && !((step-1)%freq) )
  {
    BigReal scalarPressure = controlPressure * Vector(1,1,1) / 3.;
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
    Vector &strainRate = langevinPiston_strainRate;
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
	strainRate += f2 * random->gaussian_vector();
      else
	strainRate += f2 * random->gaussian() * Vector(1,1,1);
#ifdef DEBUG_PRESSURE
      iout << iINFO << "applying langevin, strain rate: " << strainRate << "\n";
#endif
    }

    strainRate += ( 0.5 * dt * cellDims * state->lattice.volume() / mass ) *
	( controlPressure - Vector(1,1,1)*simParams->langevinPistonTarget );

#ifdef DEBUG_PRESSURE
    iout << iINFO << "integrating half step, strain rate: " << strainRate << "\n";
#endif

    if ( ! ( (step-1-slowFreq/2) % slowFreq ) )
    {
      Tensor factor;
      factor.xx = exp( dt_long * strainRate.x );
      factor.yy = exp( dt_long * strainRate.y );
      factor.zz = exp( dt_long * strainRate.z );
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
    Vector &strainRate = langevinPiston_strainRate;
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
	( controlPressure - Vector(1,1,1)*simParams->langevinPistonTarget );

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
	strainRate += f2 * random->gaussian_vector();
      else
	strainRate += f2 * random->gaussian() * Vector(1,1,1);
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
    Vector virial;
    Vector virial_normal;
    Vector virial_nbond;
    Vector virial_slow;
    Vector altVirial_normal;
    Vector altVirial_nbond;
    Vector altVirial_slow;
    Vector intVirial;
    Vector intVirial_normal;
    Vector intVirial_nbond;
    Vector intVirial_slow;
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

    virial_normal.x = reduction->item(REDUCTION_VIRIAL_NORMAL_X);
    virial_normal.y = reduction->item(REDUCTION_VIRIAL_NORMAL_Y);
    virial_normal.z = reduction->item(REDUCTION_VIRIAL_NORMAL_Z);
    virial_nbond.x = reduction->item(REDUCTION_VIRIAL_NBOND_X);
    virial_nbond.y = reduction->item(REDUCTION_VIRIAL_NBOND_Y);
    virial_nbond.z = reduction->item(REDUCTION_VIRIAL_NBOND_Z);
    virial_slow.x = reduction->item(REDUCTION_VIRIAL_SLOW_X);
    virial_slow.y = reduction->item(REDUCTION_VIRIAL_SLOW_Y);
    virial_slow.z = reduction->item(REDUCTION_VIRIAL_SLOW_Z);

    altVirial_normal.x = reduction->item(REDUCTION_ALT_VIRIAL_NORMAL_X);
    altVirial_normal.y = reduction->item(REDUCTION_ALT_VIRIAL_NORMAL_Y);
    altVirial_normal.z = reduction->item(REDUCTION_ALT_VIRIAL_NORMAL_Z);
    altVirial_nbond.x = reduction->item(REDUCTION_ALT_VIRIAL_NBOND_X);
    altVirial_nbond.y = reduction->item(REDUCTION_ALT_VIRIAL_NBOND_Y);
    altVirial_nbond.z = reduction->item(REDUCTION_ALT_VIRIAL_NBOND_Z);
    altVirial_slow.x = reduction->item(REDUCTION_ALT_VIRIAL_SLOW_X);
    altVirial_slow.y = reduction->item(REDUCTION_ALT_VIRIAL_SLOW_Y);
    altVirial_slow.z = reduction->item(REDUCTION_ALT_VIRIAL_SLOW_Z);

    intVirial_normal.x = reduction->item(REDUCTION_INT_VIRIAL_NORMAL_X);
    intVirial_normal.y = reduction->item(REDUCTION_INT_VIRIAL_NORMAL_Y);
    intVirial_normal.z = reduction->item(REDUCTION_INT_VIRIAL_NORMAL_Z);
    intVirial_nbond.x = reduction->item(REDUCTION_INT_VIRIAL_NBOND_X);
    intVirial_nbond.y = reduction->item(REDUCTION_INT_VIRIAL_NBOND_Y);
    intVirial_nbond.z = reduction->item(REDUCTION_INT_VIRIAL_NBOND_Z);
    intVirial_slow.x = reduction->item(REDUCTION_INT_VIRIAL_SLOW_X);
    intVirial_slow.y = reduction->item(REDUCTION_INT_VIRIAL_SLOW_Y);
    intVirial_slow.z = reduction->item(REDUCTION_INT_VIRIAL_SLOW_Z);

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
      pressure = Vector(0.,0.,0.);
      groupPressure = Vector(0.,0.,0.);
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
      controlPressure_normal = ( controlPressure_normal *
	Vector(1,1,1) / 3. ) * Vector(1,1,1);
      controlPressure_nbond = ( controlPressure_nbond *
	Vector(1,1,1) / 3. ) * Vector(1,1,1);
      controlPressure_slow = ( controlPressure_slow *
	Vector(1,1,1) / 3. ) * Vector(1,1,1);
      controlPressure = ( controlPressure *
	Vector(1,1,1) / 3. ) * Vector(1,1,1);
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
      NAMD_die("BUG ALERT: Bad global atom count!\n");

    checksum = reduction->item(REDUCTION_COMPUTE_CHECKSUM);
    if ( ((int)checksum) != computeChecksum ) {
      if ( computeChecksum )
        NAMD_die("BUG ALERT: Bad global compute count!\n");
      else
        computeChecksum = ((int)checksum);
    }

    checksum = reduction->item(REDUCTION_BOND_CHECKSUM);
    if ( ((int)checksum) != molecule->numCalcBonds )
      NAMD_die("BUG ALERT: Bad global bond count!\n");

    checksum = reduction->item(REDUCTION_ANGLE_CHECKSUM);
    if ( ((int)checksum) != molecule->numCalcAngles )
      NAMD_die("BUG ALERT: Bad global angle count!\n");

    checksum = reduction->item(REDUCTION_DIHEDRAL_CHECKSUM);
    if ( ((int)checksum) != molecule->numCalcDihedrals )
      NAMD_die("BUG ALERT: Bad global dihedral count!\n");

    checksum = reduction->item(REDUCTION_IMPROPER_CHECKSUM);
    if ( ((int)checksum) != molecule->numCalcImpropers )
      NAMD_die("BUG ALERT: Bad global improper count!\n");

    checksum = reduction->item(REDUCTION_EXCLUSION_CHECKSUM);
    if ( ((int)checksum) != molecule->numCalcExclusions )
      NAMD_die("BUG ALERT: Bad global exclusion count!\n");

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
  
    // Write out eXtended System Trajectory (XST) file
    if ( simParameters->xstFrequency != -1 && step == simParams->firstTimestep )
    {
      xstFile.open(simParameters->xstFilename);
      xstFile << "# NAMD extended system trajectory file" << endl;
      xstFile << "#$LABELS step a_x b_y c_z o_x o_y o_z";
      if ( simParameters->langevinPistonOn ) {
        xstFile << " s_x s_y s_z";
      }
      xstFile << endl;
    }
    if ( simParameters->xstFrequency != -1 &&
         ! ( step % simParameters->xstFrequency ) )
    {
      xstFile << step;
      xstFile << " " << lattice.a() << " " << lattice.b() << " " << lattice.c() << " " << lattice.origin().x << " " << lattice.origin().y << " " << lattice.origin().z;
      if ( simParameters->langevinPistonOn ) {
	xstFile << " " << langevinPiston_strainRate.x;
	xstFile << " " << langevinPiston_strainRate.y;
	xstFile << " " << langevinPiston_strainRate.z;
      }
      xstFile << endl;
      xstFile.flush();
    }
    if ( simParameters->xstFrequency != -1 && step == simParams->N )
    {
      xstFile.close();
    }

    // Write out eXtended System Configuration (XSC) files
    //  Output a restart file
    if ( (simParams->restartFrequency != -1) &&
         ((step % simParams->restartFrequency) == 0) &&
         (step != simParams->firstTimestep) )
    {
      char fname[140];
      strcpy(fname, simParams->restartFilename);
      strcat(fname, ".xsc");
      char bfname[140];
      strcpy(bfname, simParams->restartFilename);
      strcat(bfname, ".xsc.BAK");
      rename(fname,bfname);
      ofstream xscFile(fname);
      xscFile << "# NAMD extended system configuration file" << endl;
      xscFile << "#$LABELS step a_x b_y c_z o_x o_y o_z";
      if ( simParameters->langevinPistonOn ) {
	xscFile << " s_x s_y s_z";
      }
      xscFile << endl;
      xscFile << step;
      xscFile << " " << lattice.a() << " " << lattice.b() << " " << lattice.c() << " " << lattice.origin().x << " " << lattice.origin().y << " " << lattice.origin().z;
      if ( simParameters->langevinPistonOn ) {
	xscFile << " " << langevinPiston_strainRate.x;
	xscFile << " " << langevinPiston_strainRate.y;
	xscFile << " " << langevinPiston_strainRate.z;
      }
      xscFile << endl;
    }
    //  Output final coordinates
    if (step == simParams->N)
    {
      static char fname[140];
      strcpy(fname, simParams->outputFilename);
      strcat(fname, ".xsc");
      ofstream xscFile(fname);
      xscFile << "# NAMD extended system configuration file" << endl;
      xscFile << "#$LABELS step a_x b_y c_z o_x o_y o_z";
      if ( simParameters->langevinPistonOn ) {
	xscFile << " s_x s_y s_z";
      }
      xscFile << endl;
      xscFile << step;
      xscFile << " " << lattice.a() << " " << lattice.b() << " " << lattice.c() << " " << lattice.origin().x << " " << lattice.origin().y << " " << lattice.origin().z;
      if ( simParameters->langevinPistonOn ) {
	xscFile << " " << langevinPiston_strainRate.x;
	xscFile << " " << langevinPiston_strainRate.y;
	xscFile << " " << langevinPiston_strainRate.z;
      }
      xscFile << endl;
    }

    if (step == 32)
	startBenchTime = CmiWallTimer();
    else if (step == 64)
	iout << iINFO << "Benchmark time per step: "
	     << (CmiWallTimer() - startBenchTime) / 32. << "\n" << endi;

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

#ifdef NAMDCCS
     char webout[80];
     sprintf(webout,"%d %d %d %d",(int)totalEnergy,
	     (int)(totalEnergy - kineticEnergy),
	     (int)kineticEnergy,(int)temperature);
     CApplicationDepositNode0Data(webout);
     CkPrintf("Depositing %s\n",webout);
#endif

    if ( volume != 0. )
    {
	if ( printAtomicPressure ) {
	  iout << FORMAT(pressure*Vector(1,1,1)*PRESSUREFACTOR/3.);
	}
	iout << FORMAT(groupPressure*Vector(1,1,1)*PRESSUREFACTOR/3.);
	iout << FORMAT(volume);
    }

    if (simParams->SMDOn) {
      iout << FORMAT(smdEnergy);
    }

    iout << "\n" << endi;
}

void Controller::enqueueCollections(int timestep)
{
  if ( Output::coordinateNeeded(timestep) )
    collection->enqueuePositions(timestep);
  if ( Output::velocityNeeded(timestep) )
    collection->enqueueVelocities(timestep);
}

void Controller::terminate(void) {
  Node::Object()->enableHaltBarrier();
  CthFree(thread);
  CthSuspend();
}

