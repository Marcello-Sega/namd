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
    broadcast->positionRescaleFactor.publish(step,Vector(1,1,1)*factor);
    state->lattice.rescale(Vector(1,1,1)*factor);
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
      Vector factor;
      factor.x = exp( dt_long * strainRate.x );
      factor.y = exp( dt_long * strainRate.y );
      factor.z = exp( dt_long * strainRate.z );
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
	 electEnergy + ljEnergy + kineticEnergy + boundaryEnergy +
	 miscEnergy + smdEnergy;

    if ( node->simParameters->outputMomenta &&
         ! ( step % node->simParameters->outputMomenta ) )
    {
      iout << "MOMENTUM: " << step 
           << " P: " << momentum
           << " L: " << angularMomentum
           << "\n" << endi;
    }

#ifdef MDCOMM
    if ( node->simParameters->vmdFrequency != -1 &&
         ! ( step % node->simParameters->vmdFrequency ) )
    {
      BigReal energies[8];
      energies[0] = bondEnergy;
      energies[1] = angleEnergy;
      energies[2] = dihedralEnergy;
      energies[3] = improperEnergy;
      energies[4] = electEnergy;
      energies[5] = ljEnergy;
      energies[7] = kineticEnergy;
      Node::Object()->output->
	gather_vmd_energies(step,energies,temperature,totalEnergy);
    }
#endif  // MDCOMM

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
    if ( simParams->rigidBonds != RIGID_NONE ) { printAtomicPressure = 0; }
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
    iout << FORMAT(electEnergy);
    iout << FORMAT(ljEnergy);
    iout << FORMAT(boundaryEnergy);
    iout << FORMAT(miscEnergy);
    iout << FORMAT(kineticEnergy);
    iout << FORMAT(totalEnergy);
    iout << FORMAT(temperature);

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


/***************************************************************************
 * RCS INFORMATION:
 *
 *	$RCSfile $
 *	$Author $	$Locker:  $		$State: Exp $
 *	$Revision: 1.1075 $	$Date: 1999/07/22 15:39:41 $
 *
 ***************************************************************************
 * REVISION HISTORY:
 *
 * $Log: Controller.C,v $
 * Revision 1.1075  1999/07/22 15:39:41  jim
 * Eliminated last remnants of non-reentrant rand48 calls.
 *
 * Revision 1.1074  1999/07/08 21:26:39  jim
 * Eliminated compiler warnings.
 *
 * Revision 1.1073  1999/07/06 20:32:42  jim
 * Eliminated warnings from new generation of picky compilers.
 *
 * Revision 1.1072  1999/06/21 16:15:31  jim
 * Improved scripting, run now ends and generates output.
 *
 * Revision 1.1071  1999/06/17 17:05:44  jim
 * Renamed seq to step in most places.  Now has meaning only to user.
 *
 * Revision 1.1070  1999/06/17 15:46:14  jim
 * Completely rewrote reduction system to eliminate need for sequence numbers.
 *
 * Revision 1.1069  1999/06/03 16:50:08  jim
 * Added simplified interface to ComputeGlobal mechanism.
 *
 * Revision 1.1068  1999/06/02 15:14:19  jim
 * Now waits for output files to be written before halting.
 *
 * Revision 1.1067  1999/05/26 22:23:53  jim
 * Added basic Tcl scripting, fixed bugs in broadcasts.
 *
 * Revision 1.1066  1999/05/11 23:56:29  brunner
 * Changes for new charm version
 *
 * Revision 1.1065  1999/03/22 05:43:58  jim
 * Fixed bug in langevinPiston routine.
 *
 * Revision 1.1064  1999/03/19 23:03:00  jim
 * Fixed bugs in constant pressure code.
 *
 * Revision 1.1063  1999/03/18 02:59:11  jim
 * Eliminated pressure debugging outputs.
 *
 * Revision 1.1062  1999/03/17 21:26:32  jim
 * Switching internal nomenclature from fmaFrequency to fullElectFrequency.
 *
 * Revision 1.1061  1999/03/17 19:57:57  jim
 * Fixed logic bug in reassignment.
 *
 * Revision 1.1060  1999/03/17 17:59:23  jim
 * Eliminated compiler warnings and errors.
 *
 * Revision 1.1059  1999/03/12 02:08:35  jim
 * Fixed bug detection to deal with fixed atom optimizations.
 *
 * Revision 1.1058  1999/03/10 21:48:26  jim
 * Added remaining time to timing output.
 *
 * Revision 1.1057  1999/03/10 05:11:33  jim
 * Added reassignHold parameter.
 *
 * Revision 1.1056  1999/02/17 20:00:43  jim
 * Added error checking to keep Berendsen pressure from blowing up.
 *
 * Revision 1.1055  1999/01/06 22:50:30  jim
 * Anisotropic (flexible cell) Langevin Piston pressure control finished.
 *
 * Revision 1.1054  1999/01/06 19:19:19  jim
 * Broadcast and Sequencers understand anisotropic volume rescaling factors.
 *
 * Revision 1.1053  1999/01/06 00:56:23  jim
 * All compute objects except DPMTA now return diagonal of virial tensor.
 *
 * Revision 1.1052  1998/12/30 21:49:07  jim
 * Fixed bugs in extended system file output.
 *
 * Revision 1.1051  1998/12/07 03:54:29  jim
 * Constant pressure should work with multiple timestepping.
 * Still needs some testing.  Some debug code still enabled.
 *
 * Revision 1.1050  1998/11/30 04:10:25  krishnan
 * Added code to trigger the reduction manager on every timestep, if numNodes > nPatches
 *
 * Revision 1.1049  1998/11/30 03:19:14  jim
 * Fixed startup bug in algorithm.
 *
 * Revision 1.1048  1998/11/29 22:00:58  jim
 * Added group-based pressure control to work with rigidBonds.
 * New option useGroupPressure, turned on automatically if needed.
 *
 * Revision 1.1047  1998/11/18 21:17:41  jim
 * Added checksum to make sure compute objects don't go missing.
 *
 * Revision 1.1046  1998/11/01 23:25:47  jim
 * Added basic correctness checking: atom counts, etc.
 *
 * Revision 1.1045  1998/10/24 19:57:30  jim
 * Eliminated warnings generated by g++ -Wall.
 *
 * Revision 1.1044  1998/09/14 20:05:32  jim
 * Eliminated unnecessary timer calls.
 *
 * Revision 1.1043  1998/09/14 16:51:05  jim
 * Fixed timing printout.
 *
 * Revision 1.1042  1998/09/13 21:06:07  jim
 * Cleaned up output, defaults, etc.
 *
 * Revision 1.1041  1998/08/18 23:27:43  jim
 * First implementation of constant pressure.
 * Isotropic only, incompatible with multiple timestepping or SHAKE.
 *
 * Revision 1.1040  1998/08/04 04:07:21  jim
 * Added extended system file support and fixed lack of endi in SimParameters.
 *
 * Revision 1.1039  1998/08/03 15:31:19  jim
 * Added temperature reassignment.
 *
 * Revision 1.1038  1998/08/02 21:26:39  jim
 * Altered velocity rescaling to use averaged temperature.
 *
 * Revision 1.1037  1998/07/08 20:17:21  brunner
 * Initialized timers
 *
 * Revision 1.1036  1998/07/06 19:16:59  brunner
 * Changed path info in Makearch.T3E, changed patch partition equation,
 * and added timing prints to Controller
 *
 * Revision 1.1035  1998/06/18 15:25:22  jim
 * Workaround for aCC - put FORMAT()'s on separate lines.
 *
 * Revision 1.1034  1998/06/18 14:48:03  jim
 * Split virial into NORMAL, NBOND, and SLOW parts to match force classes.
 *
 * Revision 1.1033  1998/05/15 16:19:03  jim
 * Made Controller suspend during load balancing (for reduction system).
 *
 * Revision 1.1032  1998/04/14 03:19:20  jim
 * Fixed up MDCOMM code.
 *
 * Revision 1.1031  1998/04/06 16:34:07  jim
 * Added DPME (single processor only), test mode, and momenta printing.
 *
 * Revision 1.1030  1998/03/06 20:55:25  jim
 * Added temperature coupling.
 *
 * Revision 1.1029  1998/02/18 05:38:29  jim
 * RigidBonds mainly finished.  Now temperature is correct and a form
 * of Langevin dynamics works with constraints.
 *
 * Revision 1.1028  1998/01/05 20:25:29  sergei
 * added reduction->(un)subscribe(REDUCTION_SMD_ENERGY) to (con/de)structor
 * added SMD fields in printEnergies().
 *
 * Revision 1.1027  1997/09/25 22:52:45  brunner
 * I put in a 2-stage load balancing, so first Alg7 is done, then RefineOnly.
 * I also temporarily made the prints occur at every energy output.
 *
 * Revision 1.1026  1997/09/19 09:39:04  jim
 * Small tweaks for fixed atoms.
 *
 * Revision 1.1025  1997/08/26 16:40:22  brunner
 * Added a message to the CYCLE_BARRIER sync
 *
 * Revision 1.1024  1997/08/22 19:54:42  milind
 * Added user event EndOfTimeStep.
 *
 * Revision 1.1023  1997/08/22 19:27:36  brunner
 * Added cycle barrier, enabled by compiling with -DCYCLE_BARRIER
 *
 * Revision 1.1022  1997/08/13 14:52:20  milind
 * Made two #defines as inlined functions to fix a bug on solaris.
 *
 * Revision 1.1021  1997/04/21 00:18:53  jim
 * Fixed slowing-down problem caused by failing to require ALT_VIRIAL from
 * reduction system in controller when not printing energies.
 *
 * Revision 1.1020  1997/04/16 22:12:16  brunner
 * Fixed an LdbCoordinator bug, and cleaned up timing and Ldb output some.
 *
 * Revision 1.1019  1997/04/10 22:29:09  jim
 * First steps towards combining atom migration messages.
 *
 * Revision 1.1018  1997/04/10 18:44:33  nealk
 * 1. changed endl to endi on Controller.C
 * 2. identified popen() bug under HP-UX 9.  popen() occasionally (1/3 of the
 * time) includes garbage characters at the head of the stream, and may
 * close the stream prematurely.  I corrected the garbage bug.  Still need
 * to correct for the closing bug.  Ugh.
 *
 * Revision 1.1017  1997/04/10 17:28:43  brunner
 * Made time output in energy actually work.  The one in output doesn't do
 * anything.
 *
 * Revision 1.1016  1997/04/10 09:13:56  ari
 * Final debugging for compute migration / proxy creation for load balancing.
 * Lots of debug code added, mostly turned off now.
 * Fixed bug in PositionBox when Patch had no dependencies.
 * Eliminated use of cout and misuse of iout in numerous places.
 *                                            Ari & Jim
 *
 * Revision 1.1015  1997/04/08 21:08:40  jim
 * Contant pressure now correct on multiple nodes, should work with MTS.
 *
 * Revision 1.1014  1997/03/27 03:16:53  jim
 * Added code to check virial calculation, fixed problems with DPMTA and PBC's.
 *
 * Revision 1.1013  1997/03/25 16:57:49  nealk
 * Added PBC scaling to DPMTA.
 * Turned off debugging code in Controller.C.
 *
 * Revision 1.1012  1997/03/21 23:05:33  jim
 * Added Berendsen's pressure coupling method, won't work with MTS yet.
 *
 * Revision 1.1011  1997/03/19 22:44:21  jim
 * Revamped Controller/Sequencer, added velocity rescaling.
 *
 * Revision 1.1010  1997/03/19 11:54:13  ari
 * Add Broadcast mechanism.
 * Fixed RCS Log entries on files that did not have Log entries.
 * Added some register variables to Molecule and ComputeNonbondedExcl.C
 *
 *
 ***************************************************************************/
