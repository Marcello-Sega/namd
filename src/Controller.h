//-*-c++-*-
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

#ifndef CONTROLLER_H
#define CONTROLLER_H

#include "converse.h"
#include "Node.h"
#include "common.h"
#include <fstream.h>

class ControllerBroadcasts;
class NamdState;
class SimParameters;
class ReductionMgr;
class CollectionMaster;

class Controller
{
public:
    Controller(NamdState *s);
    ~Controller(void);
    void run(int numberOfCycles);             // spawn thread, etc.
    void awaken(void) { CthAwaken(thread); };

protected:
    virtual void algorithm(void);	// subclasses redefine this method

    void receivePressure(int);
    void printEnergies(int);
      int numDegFreedom;
      BigReal kineticEnergy;
      BigReal temperature;
      BigReal pressure;
      BigReal groupPressure;
    void enqueueCollections(int);
    void rescaleVelocities(int);
      BigReal rescaleVelocities_sumTemps;
      int rescaleVelocities_numTemps;
    void reassignVelocities(int);
    void tcoupleVelocities(int);
    void berendsenPressure(int);
    void langevinPiston1(int);
    void langevinPiston2(int);
      BigReal langevinPiston_strainRate;

    // void suspend(void) { CthSuspend(); };
    void terminate(void) {
	CPrintf("Controller terminating\n");
	Node::messageHomeDone();
	CthFree(thread); CthSuspend(); 
    };

    SimParameters *const simParams;	// for convenience
    int numberOfCycles;			// stores argument to run()
    NamdState *const state;		// access data in state
    ReductionMgr *const reduction;
    CollectionMaster *const collection;
    ControllerBroadcasts * broadcast;
    ofstream xstFile;

private:
    CthThread thread;
    static void threadRun(Controller*);

    double startCTime;
    double startWTime;
    double startBenchTime;

};

#endif // CONTROLLER_H


/***************************************************************************
 * RCS INFORMATION:
 *
 *	$RCSfile $
 *	$Author $	$Locker:  $		$State: Exp $
 *	$Revision: 1.1014 $	$Date: 1998/08/18 23:27:44 $
 *
 ***************************************************************************
 * REVISION HISTORY:
 *
 * $Log: Controller.h,v $
 * Revision 1.1014  1998/08/18 23:27:44  jim
 * First implementation of constant pressure.
 * Isotropic only, incompatible with multiple timestepping or SHAKE.
 *
 * Revision 1.1013  1998/08/04 04:07:21  jim
 * Added extended system file support and fixed lack of endi in SimParameters.
 *
 * Revision 1.1012  1998/08/03 15:31:19  jim
 * Added temperature reassignment.
 *
 * Revision 1.1011  1998/08/02 21:26:39  jim
 * Altered velocity rescaling to use averaged temperature.
 *
 * Revision 1.1010  1998/07/06 19:17:01  brunner
 * Changed path info in Makearch.T3E, changed patch partition equation,
 * and added timing prints to Controller
 *
 * Revision 1.1009  1998/03/31 04:55:44  jim
 * Added test mode, fixed errors in virial with full electrostatics.
 *
 * Revision 1.1008  1998/03/06 20:55:25  jim
 * Added temperature coupling.
 *
 * Revision 1.1007  1997/03/21 23:05:35  jim
 * Added Berendsen's pressure coupling method, won't work with MTS yet.
 *
 * Revision 1.1006  1997/03/19 22:44:23  jim
 * Revamped Controller/Sequencer, added velocity rescaling.
 *
 * Revision 1.1005  1997/03/19 11:54:14  ari
 * Add Broadcast mechanism.
 * Fixed RCS Log entries on files that did not have Log entries.
 * Added some register variables to Molecule and ComputeNonbondedExcl.C
 *
 ***************************************************************************/
