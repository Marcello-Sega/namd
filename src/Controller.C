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
#include "SimParameters.h"
#include "Controller.h"
#include "ReductionMgr.h"

#define MIN_DEBUG_LEVEL 4
#define DEBUGM
#include "Debug.h"

Controller::Controller(NamdState *s) :
	state(s),
	simParams(Node::Object()->simParameters),
	reduction(ReductionMgr::Object())
{
    reduction->subscribe(REDUCTION_BOND_ENERGY);
    reduction->subscribe(REDUCTION_ANGLE_ENERGY);
    reduction->subscribe(REDUCTION_DIHEDRAL_ENERGY);
    reduction->subscribe(REDUCTION_IMPROPER_ENERGY);
    reduction->subscribe(REDUCTION_ELECT_ENERGY);
    reduction->subscribe(REDUCTION_LJ_ENERGY);
    reduction->subscribe(REDUCTION_KINETIC_ENERGY);
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
}

void Controller::threadRun(Controller* arg)
{
    arg->algorithm();
}

void Controller::run(int numberOfCycles)
{
    stepsPerCycle = simParams->stepsPerCycle;
    if ( numberOfCycles ) this->numberOfCycles = numberOfCycles;
    else this->numberOfCycles = simParams->N; // / stepsPerCycle;
    thread = CthCreate((CthVoidFn)&(threadRun),(void*)(this),0);
    CthSetStrategyDefault(thread);
    CthAwaken(thread);
}

void Controller::algorithm(void)
{
    DebugM(4, "Controller algorithm active.\n");
    const int numberOfCycles = this->numberOfCycles;
    const int stepsPerCycle = this->stepsPerCycle;
    const BigReal timestep = simParams->dt;
    int step, cycle;
    int seq = 0;
    BigReal bondEnergy;
    BigReal angleEnergy;
    BigReal dihedralEnergy;
    BigReal improperEnergy;
    BigReal electEnergy;
    BigReal ljEnergy;
    BigReal kineticEnergy;
    BigReal totalEnergy;
iout << "Starting...\n" << endi;
    reduction->require(seq, REDUCTION_BOND_ENERGY, bondEnergy);
    reduction->require(seq, REDUCTION_ANGLE_ENERGY, angleEnergy);
    reduction->require(seq, REDUCTION_DIHEDRAL_ENERGY, dihedralEnergy);
    reduction->require(seq, REDUCTION_IMPROPER_ENERGY, improperEnergy);
    reduction->require(seq, REDUCTION_ELECT_ENERGY, electEnergy);
    reduction->require(seq, REDUCTION_LJ_ENERGY, ljEnergy);
    reduction->require(seq, REDUCTION_KINETIC_ENERGY, kineticEnergy);
iout << "Got it...\n" << endi;
    totalEnergy = bondEnergy + angleEnergy + dihedralEnergy + improperEnergy +
	 electEnergy + ljEnergy + kineticEnergy;
    iout << "ENERGY[" << seq << "] = { " <<
	"bond: " << bondEnergy << ", " << 
	"angle: " << angleEnergy << ", " << 
	"dihedral: " << dihedralEnergy << ", " << 
	"improper: " << improperEnergy << ", " << 
	"elect: " << electEnergy << ", " << 
	"lj: " << ljEnergy << ", " << 
	"kinetic: " << kineticEnergy << ", " << 
	"total: " << totalEnergy << " }\n" << endi;
    ++seq;
    for ( cycle = 0; cycle < numberOfCycles; ++cycle )
    {
        // for ( step = 0; step < stepsPerCycle; ++step )
        // {
    reduction->require(seq, REDUCTION_BOND_ENERGY, bondEnergy);
    reduction->require(seq, REDUCTION_ANGLE_ENERGY, angleEnergy);
    reduction->require(seq, REDUCTION_DIHEDRAL_ENERGY, dihedralEnergy);
    reduction->require(seq, REDUCTION_IMPROPER_ENERGY, improperEnergy);
    reduction->require(seq, REDUCTION_ELECT_ENERGY, electEnergy);
    reduction->require(seq, REDUCTION_LJ_ENERGY, ljEnergy);
    reduction->require(seq, REDUCTION_KINETIC_ENERGY, kineticEnergy);
    totalEnergy = bondEnergy + angleEnergy + dihedralEnergy + improperEnergy +
 	 electEnergy + ljEnergy + kineticEnergy;
    iout << "ENERGY[" << seq << "] = { " <<
	"bond: " << bondEnergy << ", " << 
	"angle: " << angleEnergy << ", " << 
	"dihedral: " << dihedralEnergy << ", " << 
	"improper: " << improperEnergy << ", " << 
	"elect: " << electEnergy << ", " << 
	"lj: " << ljEnergy << ", " << 
	"kinetic: " << kineticEnergy << ", " << 
	"total: " << totalEnergy << " }\n" << endi;
    ++seq;
        // }
    }
    DebugM(4, "Controller: Exiting.\n");
    terminate();
}
