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

#define MIN_DEBUG_LEVEL 3
//#define DEBUGM
#include "Debug.h"

Controller::Controller(NamdState *s) :
	state(s),
	simParams(Node::Object()->simParameters),
	reduction(ReductionMgr::Object()),
	collection(CollectionMaster::Object())
{
    if ( simParams->rescaleFreq > 0 )
    {
	velocityRescaleFactor = new
		SimpleBroadcastObject<BigReal>(velocityRescaleFactorTag);
    }
    else velocityRescaleFactor = 0;

    reduction->subscribe(REDUCTION_BOND_ENERGY);
    reduction->subscribe(REDUCTION_ANGLE_ENERGY);
    reduction->subscribe(REDUCTION_DIHEDRAL_ENERGY);
    reduction->subscribe(REDUCTION_IMPROPER_ENERGY);
    reduction->subscribe(REDUCTION_ELECT_ENERGY);
    reduction->subscribe(REDUCTION_LJ_ENERGY);
    reduction->subscribe(REDUCTION_KINETIC_ENERGY);
    reduction->subscribe(REDUCTION_BC_ENERGY);
    reduction->subscribe(REDUCTION_VIRIAL);
}

Controller::~Controller(void)
{
    delete velocityRescaleFactor;

    reduction->unsubscribe(REDUCTION_BOND_ENERGY);
    reduction->unsubscribe(REDUCTION_ANGLE_ENERGY);
    reduction->unsubscribe(REDUCTION_DIHEDRAL_ENERGY);
    reduction->unsubscribe(REDUCTION_IMPROPER_ENERGY);
    reduction->unsubscribe(REDUCTION_ELECT_ENERGY);
    reduction->unsubscribe(REDUCTION_LJ_ENERGY);
    reduction->unsubscribe(REDUCTION_KINETIC_ENERGY);
    reduction->unsubscribe(REDUCTION_BC_ENERGY);
    reduction->unsubscribe(REDUCTION_VIRIAL);
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
    thread = CthCreate((CthVoidFn)&(threadRun),(void*)(this),0);
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
    }

    terminate();
}


void Controller::rescaleVelocities(int step)
{
  const int rescaleFreq = simParams->rescaleFreq;
  if ( rescaleFreq > 0 && !(step%rescaleFreq) )
  {
    const BigReal rescaleTemp = simParams->rescaleTemp;
    BigReal factor = sqrt(rescaleTemp/temperature);
    velocityRescaleFactor->publish(step,factor);
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

    int numDegFreedom = 3 * node->molecule->numAtoms;
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

    reduction->require(seq, REDUCTION_BOND_ENERGY, bondEnergy);
    reduction->require(seq, REDUCTION_ANGLE_ENERGY, angleEnergy);
    reduction->require(seq, REDUCTION_DIHEDRAL_ENERGY, dihedralEnergy);
    reduction->require(seq, REDUCTION_IMPROPER_ENERGY, improperEnergy);
    reduction->require(seq, REDUCTION_ELECT_ENERGY, electEnergy);
    reduction->require(seq, REDUCTION_LJ_ENERGY, ljEnergy);
    reduction->require(seq, REDUCTION_KINETIC_ENERGY, kineticEnergy);
    reduction->require(seq, REDUCTION_BC_ENERGY, boundaryEnergy);
    reduction->require(seq, REDUCTION_VIRIAL, virial);

    temperature = 2.0 * kineticEnergy / ( numDegFreedom * BOLTZMAN );

    totalEnergy = bondEnergy + angleEnergy + dihedralEnergy + improperEnergy +
	 electEnergy + ljEnergy + kineticEnergy + boundaryEnergy;

    if ( seq % node->simParameters->outputEnergies ) return;

    if ( (seq % (10 * node->simParameters->outputEnergies) ) == 0 )
    {
	iout << "ETITLE:     TS    BOND        ANGLE       "
	     << "DIHED       IMPRP       ELECT       VDW       "
	     << "BOUNDARY    KINETIC        TOTAL     TEMP     "
	     << "VIRIAL\n";
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
	 << FORMAT(temperature)
	 << FORMAT(virial)
	 << "\n" << endi;

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
 *	$Revision: 1.1011 $	$Date: 1997/03/19 22:44:21 $
 *
 ***************************************************************************
 * REVISION HISTORY:
 *
 * $Log: Controller.C,v $
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
