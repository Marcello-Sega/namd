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

#ifndef cbrt
  // cbrt() not in math.h on goneril
  #define cbrt(x)  pow(x,(double)(1.0/3.0))
#endif

#define MIN_DEBUG_LEVEL 3
//#define DEBUGM
#include "Debug.h"

Controller::Controller(NamdState *s) :
	state(s),
	simParams(Node::Object()->simParameters),
	reduction(ReductionMgr::Object()),
	collection(CollectionMaster::Object()),
        startWTime(0),
        startCTime(0),
        startBenchTime(0)

{
    broadcast = new ControllerBroadcasts;

    reduction->subscribe(REDUCTION_BOND_ENERGY);
    reduction->subscribe(REDUCTION_ANGLE_ENERGY);
    reduction->subscribe(REDUCTION_DIHEDRAL_ENERGY);
    reduction->subscribe(REDUCTION_IMPROPER_ENERGY);
    reduction->subscribe(REDUCTION_ELECT_ENERGY);
    reduction->subscribe(REDUCTION_LJ_ENERGY);
    reduction->subscribe(REDUCTION_KINETIC_ENERGY);
    reduction->subscribe(REDUCTION_INT_KINETIC_ENERGY);
    reduction->subscribe(REDUCTION_BC_ENERGY);
    reduction->subscribe(REDUCTION_VIRIAL_NORMAL);
    reduction->subscribe(REDUCTION_VIRIAL_NBOND);
    reduction->subscribe(REDUCTION_VIRIAL_SLOW);
    reduction->subscribe(REDUCTION_ALT_VIRIAL_NORMAL);
    reduction->subscribe(REDUCTION_ALT_VIRIAL_NBOND);
    reduction->subscribe(REDUCTION_ALT_VIRIAL_SLOW);
    reduction->subscribe(REDUCTION_INT_VIRIAL_NORMAL);
    reduction->subscribe(REDUCTION_INT_VIRIAL_NBOND);
    reduction->subscribe(REDUCTION_INT_VIRIAL_SLOW);
    reduction->subscribe(REDUCTION_SMD_ENERGY);
    reduction->subscribe(REDUCTION_MOMENTUM_X);
    reduction->subscribe(REDUCTION_MOMENTUM_Y);
    reduction->subscribe(REDUCTION_MOMENTUM_Z);
    reduction->subscribe(REDUCTION_ANGULAR_MOMENTUM_X);
    reduction->subscribe(REDUCTION_ANGULAR_MOMENTUM_Y);
    reduction->subscribe(REDUCTION_ANGULAR_MOMENTUM_Z);
}

Controller::~Controller(void)
{
    delete broadcast;

    reduction->unsubscribe(REDUCTION_BOND_ENERGY);
    reduction->unsubscribe(REDUCTION_ANGLE_ENERGY);
    reduction->unsubscribe(REDUCTION_DIHEDRAL_ENERGY);
    reduction->unsubscribe(REDUCTION_IMPROPER_ENERGY);
    reduction->unsubscribe(REDUCTION_ELECT_ENERGY);
    reduction->unsubscribe(REDUCTION_LJ_ENERGY);
    reduction->unsubscribe(REDUCTION_KINETIC_ENERGY);
    reduction->unsubscribe(REDUCTION_INT_KINETIC_ENERGY);
    reduction->unsubscribe(REDUCTION_BC_ENERGY);
    reduction->unsubscribe(REDUCTION_VIRIAL_NORMAL);
    reduction->unsubscribe(REDUCTION_VIRIAL_NBOND);
    reduction->unsubscribe(REDUCTION_VIRIAL_SLOW);
    reduction->unsubscribe(REDUCTION_ALT_VIRIAL_NORMAL);
    reduction->unsubscribe(REDUCTION_ALT_VIRIAL_NBOND);
    reduction->unsubscribe(REDUCTION_ALT_VIRIAL_SLOW);
    reduction->unsubscribe(REDUCTION_INT_VIRIAL_NORMAL);
    reduction->unsubscribe(REDUCTION_INT_VIRIAL_NBOND);
    reduction->unsubscribe(REDUCTION_INT_VIRIAL_SLOW);
    reduction->unsubscribe(REDUCTION_SMD_ENERGY);
    reduction->unsubscribe(REDUCTION_MOMENTUM_X);
    reduction->unsubscribe(REDUCTION_MOMENTUM_Y);
    reduction->unsubscribe(REDUCTION_MOMENTUM_Z);
    reduction->unsubscribe(REDUCTION_ANGULAR_MOMENTUM_X);
    reduction->unsubscribe(REDUCTION_ANGULAR_MOMENTUM_Y);
    reduction->unsubscribe(REDUCTION_ANGULAR_MOMENTUM_Z);
}

void Controller::threadRun(Controller* arg)
{
    arg->algorithm();
}

void Controller::run(int numberOfCycles)
{
    this->numberOfCycles = numberOfCycles;
    if ( numberOfCycles ) 
      NAMD_die("Sorry, Controller::run() does not support an argument.\n");

    // create a Thread and invoke it
    DebugM(4, "Starting thread in controller on this=" << this << "\n");
    thread = CthCreate((CthVoidFn)&(threadRun),(void*)(this),CTRL_STK_SZ);
    CthSetStrategyDefault(thread);
    awaken();
}

extern int eventEndOfTimeStep;
extern "C" void trace_user_event(int event);

void Controller::algorithm(void)
{
    int step = simParams->firstTimestep;

    const int numberOfSteps = simParams->N;
    const int stepsPerCycle = simParams->stepsPerCycle;
    const BigReal timestep = simParams->dt;

    for ( ; step <= numberOfSteps; ++step )
    {
        enqueueCollections(step);
        trace_user_event(eventEndOfTimeStep);
        printEnergies(step);
        rescaleVelocities(step);
	tcoupleVelocities(step);
	berendsenPressure(step);
#ifdef CYCLE_BARRIER
	if (!((step+1) % stepsPerCycle))
	{
	  broadcast->cycleBarrier.publish(step,1);
	  CPrintf("Cycle time at sync Wall: %f CPU %f\n",
		  CmiWallTimer(),CmiTimer());
	}
#endif
	if ( LdbCoordinator::Object()->balanceNow(step) ) {
	  LdbCoordinator::Object()->rebalance(this);
	}
    }

    terminate();
}

void Controller::berendsenPressure(int step)
{
  const int freq = simParams->berendsenPressureFreq;
  if ( simParams->berendsenPressureOn && !(step%freq) )
  {
    BigReal factor = pressure - simParams->berendsenPressureTarget;
    factor *= simParams->berendsenPressureCompressibility;
    factor *= ( simParams->dt * freq );
    factor /= simParams->berendsenPressureRelaxationTime;
    factor += 1.0;
    factor = cbrt(factor);
    broadcast->positionRescaleFactor.publish(step,factor);
    state->lattice.rescale(factor);
  }
}

void Controller::rescaleVelocities(int step)
{
  const int rescaleFreq = simParams->rescaleFreq;
  if ( rescaleFreq > 0 && !(step%rescaleFreq) )
  {
    const BigReal rescaleTemp = simParams->rescaleTemp;
    BigReal factor = sqrt(rescaleTemp/temperature);
    broadcast->velocityRescaleFactor.publish(step,factor);
    iout << "RESCALING VELOCITIES AT STEP " << step
         << " TO " << rescaleTemp << " KELVIN.\n" << endi;
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

void Controller::printEnergies(int seq)
{
    Node *node = Node::Object();
    Molecule *molecule = node->molecule;
    SimParameters *simParameters = node->simParameters;
    Lattice &lattice = state->lattice;

    int numAtoms = molecule->numAtoms;
    int numDegFreedom = 3 * numAtoms;
    int numFixedAtoms = molecule->numFixedAtoms;
    if ( numFixedAtoms ) numDegFreedom -= 3 * numFixedAtoms;
    if ( ! ( numFixedAtoms || molecule->numConstraints
	|| simParameters->comMove || simParameters->langevinOn ) ) {
      numDegFreedom -= 3;
    }
    int numRigidBonds = molecule->numRigidBonds;
    int numFixedRigidBonds = molecule->numFixedRigidBonds;
    numDegFreedom -= ( numRigidBonds - numFixedRigidBonds );

    BigReal bondEnergy;
    BigReal angleEnergy;
    BigReal dihedralEnergy;
    BigReal improperEnergy;
    BigReal electEnergy;
    BigReal ljEnergy;
    BigReal kineticEnergy;
    BigReal intKineticEnergy;
    BigReal boundaryEnergy;
    BigReal tmpVirial;
    BigReal virial;
    BigReal altVirial;
    BigReal intVirial;
    BigReal smdEnergy;
    BigReal totalEnergy;
    BigReal volume;
    Vector momentum;
    Vector angularMomentum;

    reduction->require(seq, REDUCTION_BOND_ENERGY, bondEnergy);
    reduction->require(seq, REDUCTION_ANGLE_ENERGY, angleEnergy);
    reduction->require(seq, REDUCTION_DIHEDRAL_ENERGY, dihedralEnergy);
    reduction->require(seq, REDUCTION_IMPROPER_ENERGY, improperEnergy);
    reduction->require(seq, REDUCTION_ELECT_ENERGY, electEnergy);
    reduction->require(seq, REDUCTION_LJ_ENERGY, ljEnergy);
    reduction->require(seq, REDUCTION_KINETIC_ENERGY, kineticEnergy);
    reduction->require(seq, REDUCTION_INT_KINETIC_ENERGY, intKineticEnergy);
    reduction->require(seq, REDUCTION_BC_ENERGY, boundaryEnergy);
    reduction->require(seq, REDUCTION_SMD_ENERGY, smdEnergy);
    virial = 0;
    reduction->require(seq, REDUCTION_VIRIAL_NORMAL, tmpVirial);
    virial += tmpVirial;
    reduction->require(seq, REDUCTION_VIRIAL_NBOND, tmpVirial);
    virial += tmpVirial;
    reduction->require(seq, REDUCTION_VIRIAL_SLOW, tmpVirial);
    virial += tmpVirial;
    virial /= 3.;  // virial submitted is wrong by factor of 3
    altVirial = 0;
    reduction->require(seq, REDUCTION_ALT_VIRIAL_NORMAL, tmpVirial);
    altVirial += tmpVirial;
    reduction->require(seq, REDUCTION_ALT_VIRIAL_NBOND, tmpVirial);
    altVirial += tmpVirial;
    reduction->require(seq, REDUCTION_ALT_VIRIAL_SLOW, tmpVirial);
    altVirial += tmpVirial;
    altVirial /= 3.;  // virial submitted is wrong by factor of 3
    intVirial = 0;
    reduction->require(seq, REDUCTION_INT_VIRIAL_NORMAL, tmpVirial);
    intVirial += tmpVirial;
    reduction->require(seq, REDUCTION_INT_VIRIAL_NBOND, tmpVirial);
    intVirial += tmpVirial;
    reduction->require(seq, REDUCTION_INT_VIRIAL_SLOW, tmpVirial);
    intVirial += tmpVirial;
    intVirial /= 3.;  // virial submitted is wrong by factor of 3

    reduction->require(seq, REDUCTION_MOMENTUM_X, momentum.x);
    reduction->require(seq, REDUCTION_MOMENTUM_Y, momentum.y);
    reduction->require(seq, REDUCTION_MOMENTUM_Z, momentum.z);
    reduction->require(seq, REDUCTION_ANGULAR_MOMENTUM_X, angularMomentum.x);
    reduction->require(seq, REDUCTION_ANGULAR_MOMENTUM_Y, angularMomentum.y);
    reduction->require(seq, REDUCTION_ANGULAR_MOMENTUM_Z, angularMomentum.z);

    temperature = 2.0 * kineticEnergy / ( numDegFreedom * BOLTZMAN );

    BigReal groupPressure;
    if ( (volume=lattice.volume()) != 0. )
    {
      pressure = ( numAtoms * BOLTZMAN * temperature + virial ) / volume;
      groupPressure = ( (2./3.)*( kineticEnergy - intKineticEnergy ) +
                        ( virial - intVirial ) ) / volume;
    }
    else
    {
      pressure = 0.;
      groupPressure = 0.;
    }

    totalEnergy = bondEnergy + angleEnergy + dihedralEnergy + improperEnergy +
	 electEnergy + ljEnergy + kineticEnergy + boundaryEnergy + smdEnergy;

    if ( node->simParameters->outputMomenta &&
         ! ( seq % node->simParameters->outputMomenta ) )
    {
      iout << "MOMENTA: " << seq 
           << " P=" << momentum
           << " L=" << angularMomentum
           << "\n" << endi;
    }

#ifdef MDCOMM
    if ( node->simParameters->vmdFrequency != -1 &&
         ! ( seq % node->simParameters->vmdFrequency ) )
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
	gather_vmd_energies(seq,energies,temperature,totalEnergy);
    }
#endif  // MDCOMM

    // NO CALCULATIONS OR REDUCTIONS BEYOND THIS POINT!!!
    if ( seq % node->simParameters->outputEnergies ) return;
    // ONLY OUTPUT SHOULD OCCUR BELOW THIS LINE!!!

    const double endWTime = CmiWallTimer();
    const double endCTime = CmiTimer();

    if (seq == 32)
	startBenchTime = endWTime;
    else if (seq == 64)
	iout << iINFO << "Benchmark time: "
	     << (endWTime - startBenchTime) / 32. << "\n" << endi;

    const double elapsedW = 
	(endWTime - startWTime) / node->simParameters->outputEnergies;
    const double elapsedC = 
	(endCTime - startCTime) / node->simParameters->outputEnergies;

    startWTime = endWTime;
    startCTime = endCTime;

    iout << iINFO
    	 << "Elapsed time CPU: " << endCTime << " Wall: " 
    	 << endWTime << "\n" << endi;

    iout << iINFO
    	 << "Time per step CPU: " << elapsedC << " Wall: " 
    	 << elapsedW << "\n" << endi;

    if ( (seq % (10 * node->simParameters->outputEnergies) ) == 0 )
    {
	iout << "ETITLE:     TS    BOND        ANGLE       "
	     << "DIHED       IMPRP       ELECT       VDW       "
	     << "BOUNDARY    KINETIC        TOTAL     TEMP";
	if ( volume != 0. ) iout << "     PRESSURE    GPRESSURE    VOLUME";
	if (simParams->SMDOn) iout << "     SMD";
	iout << "\n" << endi;
    }

    // N.B.  HP's aCC compiler merges FORMAT calls in the same expression.
    //       Need separate statements because data returned in static array.

    iout << ETITLE(seq);
    iout << FORMAT(bondEnergy);
    iout << FORMAT(angleEnergy);
    iout << FORMAT(dihedralEnergy);
    iout << FORMAT(improperEnergy);
    iout << FORMAT(electEnergy);
    iout << FORMAT(ljEnergy);
    iout << FORMAT(boundaryEnergy);
    iout << FORMAT(kineticEnergy);
    iout << FORMAT(totalEnergy);
    iout << FORMAT(temperature);

    if ( volume != 0. )
    {
	iout << FORMAT(pressure*PRESSUREFACTOR);
	iout << FORMAT(groupPressure*PRESSUREFACTOR);
	iout << FORMAT(volume);
    }

    if (simParams->SMDOn) {
      iout << FORMAT(smdEnergy);
    }

    iout << "\n" << endi;

    DebugM(4,"step: " << seq << " virial: " << virial
		<< " altVirial: " << altVirial << "\n");

}

void Controller::enqueueCollections(int timestep)
{
  if ( Output::coordinateNeeded(timestep) )
    collection->enqueuePositions(timestep);
  if ( Output::velocityNeeded(timestep) )
    collection->enqueueVelocities(timestep);
}

/***************************************************************************
 * RCS INFORMATION:
 *
 *	$RCSfile $
 *	$Author $	$Locker:  $		$State: Exp $
 *	$Revision: 1.1037 $	$Date: 1998/07/08 20:17:21 $
 *
 ***************************************************************************
 * REVISION HISTORY:
 *
 * $Log: Controller.C,v $
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
