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


/***************************************************************************
 * RCS INFORMATION:
 *
 *	$RCSfile $
 *	$Author $	$Locker:  $		$State: Exp $
 *	$Revision: 1.1 $	$Date: 1998/03/31 04:55:50 $
 *
 ***************************************************************************
 * REVISION HISTORY:
 *
 * $Log: TestSequencer.h,v $
 * Revision 1.1  1998/03/31 04:55:50  jim
 * Added test mode, fixed errors in virial with full electrostatics.
 *
 *
 ***************************************************************************/
