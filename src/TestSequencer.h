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

