/**
***  Copyright (c) 1995, 1996, 1997, 1998, 1999, 2000 by
***  The Board of Trustees of the University of Illinois.
***  All rights reserved.
**/

#ifndef COMPUTEEFIELD_H
#define COMPUTEEFIELD_H

#include "ComputePatch.h"
#include "ReductionMgr.h"

class ComputeEField : public ComputePatch
{

public:
	ComputeEField(ComputeID c, PatchID pid); 	//  Constructor
	virtual ~ComputeEField();			//  Destructor

	virtual void doForce(CompAtom* p, Results* r);

	SubmitReduction *reduction;

};

#endif







