/**
***  Copyright (c) 1995, 1996, 1997, 1998, 1999, 2000 by
***  The Board of Trustees of the University of Illinois.
***  All rights reserved.
**/

#ifndef TESTSEQUENCER_H
#define TESTSEQUENCER_H

#include "Sequencer.h"

class TestSequencer : public Sequencer
{
public:
    TestSequencer(HomePatch *p);
    ~TestSequencer(void);

protected:
    virtual void algorithm(void);	// subclasses redefine this method

    void translatePosition(BigReal,BigReal,BigReal);

};

#endif // TESTSEQUENCER_H

