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

class HomePatch;
class SimParameters;
class ReductionMgr;
class CollectionMgr;
template <class T>
class SimpleBroadcastObject;

class Sequencer
{
    friend HomePatch;
public:
    Sequencer(HomePatch *p);
    ~Sequencer(void);
    void run(int numberOfCycles);             // spawn thread, etc.
    void awaken(void) { CthAwaken(thread); }

protected:
    virtual void algorithm(void);	// subclasses redefine this method

    void submitCollections(int);
    void rescaleVelocities(int);

    void suspend(void) { CthSuspend(); }
    void terminate(void);
    SimParameters *const simParams;	// for convenience
    int numberOfCycles;			// stores argument to run()
    HomePatch *const patch;		// access methods in patch
    ReductionMgr *const reduction;
    CollectionMgr *const collection;

private:
    CthThread thread;
    static void threadRun(Sequencer*);

    SimpleBroadcastObject<BigReal> * velocityRescaleFactor;
};

#endif // SEQUENCER_H


/***************************************************************************
 * RCS INFORMATION:
 *
 *	$RCSfile $
 *	$Author $	$Locker:  $		$State: Exp $
 *	$Revision: 1.1006 $	$Date: 1997/03/19 22:44:26 $
 *
 ***************************************************************************
 * REVISION HISTORY:
 *
 * $Log: Sequencer.h,v $
 * Revision 1.1006  1997/03/19 22:44:26  jim
 * Revamped Controller/Sequencer, added velocity rescaling.
 *
 * Revision 1.1005  1997/03/19 11:54:56  ari
 * Add Broadcast mechanism.
 * Fixed RCS Log entries on files that did not have Log entries.
 * Added some register variables to Molecule and ComputeNonbondedExcl.C
 *
 ***************************************************************************/
