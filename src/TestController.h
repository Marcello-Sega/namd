/**
***  Copyright (c) 1995, 1996, 1997, 1998, 1999, 2000 by
***  The Board of Trustees of the University of Illinois.
***  All rights reserved.
**/

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

