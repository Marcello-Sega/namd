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

