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

#ifndef TESTCONTROLLER_H
#define TESTCONTROLLER_H

#include "Controller.h"

class TestController : public Controller
{
public:
    TestController(NamdState *s);
    ~TestController(void);

protected:
    virtual void algorithm(void);	// subclasses redefine this method
    void berendsenPressure(int);

};

#endif // TESTCONTROLLER_H


/***************************************************************************
 * RCS INFORMATION:
 *
 *	$RCSfile $
 *	$Author $	$Locker:  $		$State: Exp $
 *	$Revision: 1.1 $	$Date: 1998/03/31 04:55:49 $
 *
 ***************************************************************************
 * REVISION HISTORY:
 *
 * $Log: TestController.h,v $
 * Revision 1.1  1998/03/31 04:55:49  jim
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
