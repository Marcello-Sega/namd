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

class NamdState;
class SimParameters;
class ReductionMgr;

class Controller
{
public:
    Controller(NamdState *s);
    ~Controller(void);
    void run(int numberOfCycles);             // spawn thread, etc.
    void awaken(void) { CthAwaken(thread); };

protected:
    virtual void algorithm(void);	// subclasses redefine this method

    void suspend(void) { CthSuspend(); };
    void terminate(void) { CthFree(thread); CthSuspend(); };
    SimParameters *const simParams;	// for convenience
    int numberOfCycles;			// stores argument to run()
    int stepsPerCycle;			// stores info from run()
    NamdState *const state;		// access data in state
    ReductionMgr *const reduction;

private:
    CthThread thread;
    static void threadRun(Controller*);
};

#endif // SEQUENCER_H

