/**
***  Copyright (c) 1995, 1996, 1997, 1998, 1999, 2000 by
***  The Board of Trustees of the University of Illinois.
***  All rights reserved.
**/

#ifndef COMPUTEGRIDFORCE_H
#define COMPUTEGRIDFORCE_H

#include "ComputeHomePatch.h"
#include "ReductionMgr.h"
#include "GridForceGrid.h"

class ComputeGridForce : public ComputeHomePatch
{

public:
    ComputeGridForce(ComputeID c, PatchID pid); 	//  Constructor
    virtual ~ComputeGridForce();			//  Destructor
    
    void doForce(FullAtom* p, Results* r);
    
    SubmitReduction *reduction;
};

#endif
