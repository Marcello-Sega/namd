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

class Sequencer
{
public:
    Sequencer(HomePatch *p) : patch(p) { };
    ~Sequencer(void) { };
    void run(int numberOfCycles);             // spawn thread, etc.
    void awaken(void) { CthAwaken(thread); };

protected:
    void suspend(void) { CthSuspend(); };
    void terminate(void) { CthFree(thread); CthSuspend(); };
    virtual void threadRun(void);  // subclasses redefine this method
    int numberOfCycles;            // stores argument to run()
    HomePatch *const patch;        // access methods in patch

private:
    CthThread thread;
    friend void SequencerThreadRun(Sequencer*);
};

#endif // SEQUENCER_H

