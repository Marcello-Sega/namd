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
    sequence = new SimpleBroadcastObject<int>(1);

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
    stepsPerCycle = simParams->stepsPerCycle;
    if ( numberOfCycles ) 
      this->numberOfCycles = numberOfCycles;
    else 
      this->numberOfCycles = simParams->N - simParams->firstTimestep; 
      // / stepsPerCycle;
    DebugM(4, "Starting thread in controller on this=" << this << "\n");
    thread = CthCreate((CthVoidFn)&(threadRun),(void*)(this),0);
    CthSetStrategyDefault(thread);
    CthAwaken(thread);
}

void Controller::algorithm(void)
{
    const int numberOfCycles = this->numberOfCycles;
    const int stepsPerCycle = this->stepsPerCycle;
    const BigReal timestep = simParams->dt;
    const int first = simParams->firstTimestep;
    DebugM(4, "Controller algorithm active. timestep = " << timestep << "\n");
    // int step;
    int cycle;
    int seq = 0;
    enqueueCollections(seq+first);
    printEnergies(seq);
    //sequence->publish(seq,seq); // broadcast the value seq to all Sequencers
    seq++;
    for ( cycle = 0; cycle < numberOfCycles; ++cycle )
    {
        enqueueCollections(seq+first);
        sequence->publish(seq,seq);
        printEnergies(seq);
	seq++;
    }
    terminate();
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
    BigReal temperature;
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

    iout << ETITLE(seq + simParams->firstTimestep)
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
 *	$Revision: 1.1010 $	$Date: 1997/03/19 11:54:13 $
 *
 ***************************************************************************
 * REVISION HISTORY:
 *
 * $Log: Controller.C,v $
 * Revision 1.1010  1997/03/19 11:54:13  ari
 * Add Broadcast mechanism.
 * Fixed RCS Log entries on files that did not have Log entries.
 * Added some register variables to Molecule and ComputeNonbondedExcl.C
 *
 *
 ***************************************************************************/
