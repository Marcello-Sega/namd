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
    ~Sequencer(void);
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

    void rescaleVelocities(int);
    void berendsenPressure(int);
    void langevinVelocities(int);

    void terminate(void);
    SimParameters *const simParams;	// for convenience
    int numberOfCycles;			// stores argument to run()
    HomePatch *const patch;		// access methods in patch
    ReductionMgr *const reduction;
    CollectionMgr *const collection;

private:
    void rebalanceLoad(int timestep);
    CthThread thread;
    static void threadRun(Sequencer*);

    ControllerBroadcasts * broadcast;
    LdbCoordinator *ldbCoordinator;
};

#endif // SEQUENCER_H


/***************************************************************************
 * RCS INFORMATION:
 *
 *	$RCSfile $
 *	$Author $	$Locker:  $		$State: Exp $
 *	$Revision: 1.1012 $	$Date: 1998/02/17 06:39:24 $
 *
 ***************************************************************************
 * REVISION HISTORY:
 *
 * $Log: Sequencer.h,v $
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
