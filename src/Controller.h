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

template <class T>
class SimpleBroadcastObject;

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

    void printEnergies(int);
    void enqueueCollections(int);

    // void suspend(void) { CthSuspend(); };
    void terminate(void) {
	CPrintf("Controller terminating\n");
	Node::messageHomeDone();
	CthFree(thread); CthSuspend(); 
    };

    SimParameters *const simParams;	// for convenience
    int numberOfCycles;			// stores argument to run()
    int stepsPerCycle;			// stores info from run()
    NamdState *const state;		// access data in state
    ReductionMgr *const reduction;
    CollectionMaster *const collection;

private:
    CthThread thread;
    static void threadRun(Controller*);

    SimpleBroadcastObject<int> * sequence;
};

#endif // SEQUENCER_H


/***************************************************************************
 * RCS INFORMATION:
 *
 *	$RCSfile $
 *	$Author $	$Locker:  $		$State: Exp $
 *	$Revision: 1.1005 $	$Date: 1997/03/19 11:54:14 $
 *
 ***************************************************************************
 * REVISION HISTORY:
 *
 * $Log: Controller.h,v $
 * Revision 1.1005  1997/03/19 11:54:14  ari
 * Add Broadcast mechanism.
 * Fixed RCS Log entries on files that did not have Log entries.
 * Added some register variables to Molecule and ComputeNonbondedExcl.C
 *
 ***************************************************************************/
