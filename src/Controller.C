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
	collection(CollectionMaster::Object())
{
    broadcast = new ControllerBroadcasts;

    reduction->subscribe(REDUCTION_BOND_ENERGY);
    reduction->subscribe(REDUCTION_ANGLE_ENERGY);
    reduction->subscribe(REDUCTION_DIHEDRAL_ENERGY);
    reduction->subscribe(REDUCTION_IMPROPER_ENERGY);
    reduction->subscribe(REDUCTION_ELECT_ENERGY);
    reduction->subscribe(REDUCTION_LJ_ENERGY);
    reduction->subscribe(REDUCTION_KINETIC_ENERGY);
    reduction->subscribe(REDUCTION_BC_ENERGY);
    reduction->subscribe(REDUCTION_VIRIAL);
    reduction->subscribe(REDUCTION_ALT_VIRIAL);
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
    reduction->unsubscribe(REDUCTION_BC_ENERGY);
    reduction->unsubscribe(REDUCTION_VIRIAL);
    reduction->unsubscribe(REDUCTION_ALT_VIRIAL);
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

void Controller::algorithm(void)
{
    int step = simParams->firstTimestep;

    const int numberOfSteps = simParams->N;
    const int stepsPerCycle = simParams->stepsPerCycle;
    const BigReal timestep = simParams->dt;

    for ( ; step <= numberOfSteps; ++step )
    {
        enqueueCollections(step);
        printEnergies(step);
        rescaleVelocities(step);
	berendsenPressure(step);
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

void Controller::printEnergies(int seq)
{
    char tmp_string[21];
#define ETITLE(X) ( sprintf(tmp_string,"ENERGY: %6d ",X), tmp_string )
#define FORMAT(X) ( sprintf(tmp_string,"%.4f ",X), NAMD_pad(tmp_string, 12), tmp_string )

    Node *node = Node::Object();
    Lattice &lattice = state->lattice;

    int numAtoms = node->molecule->numAtoms;
    int numDegFreedom = 3 * numAtoms;
    if ( ! node->simParameters->comMove ) numDegFreedom -= 3;

    BigReal bondEnergy;
    BigReal angleEnergy;
    BigReal dihedralEnergy;
    BigReal improperEnergy;
    BigReal electEnergy;
    BigReal ljEnergy;
    BigReal kineticEnergy;
    BigReal boundaryEnergy;
    BigReal virial;
    BigReal totalEnergy;
    BigReal volume;

    reduction->require(seq, REDUCTION_BOND_ENERGY, bondEnergy);
    reduction->require(seq, REDUCTION_ANGLE_ENERGY, angleEnergy);
    reduction->require(seq, REDUCTION_DIHEDRAL_ENERGY, dihedralEnergy);
    reduction->require(seq, REDUCTION_IMPROPER_ENERGY, improperEnergy);
    reduction->require(seq, REDUCTION_ELECT_ENERGY, electEnergy);
    reduction->require(seq, REDUCTION_LJ_ENERGY, ljEnergy);
    reduction->require(seq, REDUCTION_KINETIC_ENERGY, kineticEnergy);
    reduction->require(seq, REDUCTION_BC_ENERGY, boundaryEnergy);
    reduction->require(seq, REDUCTION_VIRIAL, virial);

    virial /= 3.;  // virial submitted is wrong by factor of 3

    temperature = 2.0 * kineticEnergy / ( numDegFreedom * BOLTZMAN );

    if ( (volume=lattice.volume()) != 0. )
    {
      pressure = ( numAtoms * BOLTZMAN * temperature + virial ) / volume;
    }
    else
    {
      pressure = 0.;
    }

    totalEnergy = bondEnergy + angleEnergy + dihedralEnergy + improperEnergy +
	 electEnergy + ljEnergy + kineticEnergy + boundaryEnergy;

    if ( seq % node->simParameters->outputEnergies ) return;

    if ( (seq % (10 * node->simParameters->outputEnergies) ) == 0 )
    {
	iout << "ETITLE:     TS    BOND        ANGLE       "
	     << "DIHED       IMPRP       ELECT       VDW       "
	     << "BOUNDARY    KINETIC        TOTAL     TEMP";
	if ( volume != 0. ) iout << "     PRESSURE    VOLUME";
	iout << "\n";
    }

    iout << ETITLE(seq)
	 << FORMAT(bondEnergy)
	 << FORMAT(angleEnergy)
	 << FORMAT(dihedralEnergy)
	 << FORMAT(improperEnergy)
	 << FORMAT(electEnergy)
	 << FORMAT(ljEnergy)
	 << FORMAT(boundaryEnergy)
	 << FORMAT(kineticEnergy)
	 << FORMAT(totalEnergy)
	 << FORMAT(temperature);
    if ( volume != 0. )
    {
	iout << FORMAT(pressure*PRESSUREFACTOR)
	     << FORMAT(volume);
    }
    iout << "\n" << endi;

    BigReal altVirial;
    reduction->require(seq, REDUCTION_ALT_VIRIAL, altVirial);
    altVirial /= 3.;  // virial submitted is wrong by factor of 3
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
 *	$Revision: 1.1016 $	$Date: 1997/04/10 09:13:56 $
 *
 ***************************************************************************
 * REVISION HISTORY:
 *
 * $Log: Controller.C,v $
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
