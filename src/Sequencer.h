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

#ifndef SEQUENCER_H
#define SEQUENCER_H

#include "converse.h"
#include "PatchTypes.h"

class HomePatch;
class SimParameters;
class ReductionMgr;
class CollectionMgr;
class ControllerBroadcasts;
class LdbCoordinator;

class Sequencer
{
    friend class HomePatch;
public:
    Sequencer(HomePatch *p);
    virtual ~Sequencer(void);
    void run(int numberOfCycles);             // spawn thread, etc.
    void awaken(void) { CthAwaken(thread); }
    void suspend(void) { CthSuspend(); }

protected:
    virtual void algorithm(void);	// subclasses redefine this method

    void runComputeObjects(int migration = 0);

    void submitReductions(int);
    void submitCollections(int);

    void addForceToMomentum(BigReal, const int ftag = Results::normal);
    void addVelocityToPosition(BigReal);

    void rattle1(BigReal);
    void rattle2(BigReal,int);

    void minimizationMaximumMove(BigReal);
    void minimizationQuenchVelocity(void);

    void rescaleVelocities(int);
      int rescaleVelocities_numTemps;
    void reassignVelocities(int);
    void tcoupleVelocities(BigReal,int);
    void berendsenPressure(int);
    void langevinPiston(int);
      int slowFreq;
    void langevinVelocities(BigReal);
    void langevinVelocitiesBBK1(BigReal);
    void langevinVelocitiesBBK2(BigReal);

    void terminate(void);
    SimParameters *const simParams;	// for convenience
    int numberOfCycles;			// stores argument to run()
    HomePatch *const patch;		// access methods in patch
    ReductionMgr *const reduction;
    CollectionMgr *const collection;
    ControllerBroadcasts * broadcast;

    void rebalanceLoad(int timestep);

private:
    CthThread thread;
    static void threadRun(Sequencer*);

    LdbCoordinator *ldbCoordinator;
};

#endif // SEQUENCER_H


/***************************************************************************
 * RCS INFORMATION:
 *
 *	$RCSfile $
 *	$Author $	$Locker:  $		$State: Exp $
 *	$Revision: 1.1022 $	$Date: 1999/04/27 23:43:03 $
 *
 ***************************************************************************
 * REVISION HISTORY:
 *
 * $Log: Sequencer.h,v $
 * Revision 1.1022  1999/04/27 23:43:03  jim
 * Switched Langevin dynamics integrator to a two-part version of BBK.
 *
 * Revision 1.1021  1999/03/19 23:03:02  jim
 * Fixed bugs in constant pressure code.
 *
 * Revision 1.1020  1998/10/24 19:58:00  jim
 * Eliminated warnings generated by g++ -Wall.
 *
 * Revision 1.1019  1998/08/18 23:27:45  jim
 * First implementation of constant pressure.
 * Isotropic only, incompatible with multiple timestepping or SHAKE.
 *
 * Revision 1.1018  1998/08/03 15:31:20  jim
 * Added temperature reassignment.
 *
 * Revision 1.1017  1998/08/02 21:26:41  jim
 * Altered velocity rescaling to use averaged temperature.
 *
 * Revision 1.1016  1998/03/31 04:55:46  jim
 * Added test mode, fixed errors in virial with full electrostatics.
 *
 * Revision 1.1015  1998/03/06 20:55:26  jim
 * Added temperature coupling.
 *
 * Revision 1.1014  1998/03/06 10:25:28  jim
 * Added very basic minimizer.
 *
 * Revision 1.1013  1998/02/18 19:14:00  jim
 * Fixed Langevin dynamics, undoing changes from yesterday.
 *
 * Revision 1.1012  1998/02/17 06:39:24  jim
 * SHAKE/RATTLE (rigidBonds) appears to work!!!  Still needs langevin,
 * proper startup, and degree of freedom tracking.
 *
 * Revision 1.1011  1998/01/15 04:58:50  jim
 * Corrected "friend foo" to "friend class foo".
 *
 * Revision 1.1010  1997/03/27 20:25:52  brunner
 * Changes for LdbCoordinator, the load balance control BOC
 *
 * Revision 1.1009  1997/03/25 04:04:59  jim
 * Simplified algorithm a bit.
 *
 * Revision 1.1008  1997/03/24 01:44:02  jim
 * Added Langevin dynamics.
 *
 * Revision 1.1007  1997/03/21 23:05:43  jim
 * Added Berendsen's pressure coupling method, won't work with MTS yet.
 *
 * Revision 1.1006  1997/03/19 22:44:26  jim
 * Revamped Controller/Sequencer, added velocity rescaling.
 *
 * Revision 1.1005  1997/03/19 11:54:56  ari
 * Add Broadcast mechanism.
 * Fixed RCS Log entries on files that did not have Log entries.
 * Added some register variables to Molecule and ComputeNonbondedExcl.C
 *
 ***************************************************************************/
