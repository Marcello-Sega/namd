/**
***  Copyright (c) 1995, 1996, 1997, 1998, 1999, 2000 by
***  The Board of Trustees of the University of Illinois.
***  All rights reserved.
**/

#include "ComputeEField.h"
#include "Node.h"
#include "SimParameters.h"
#include "Patch.h"


ComputeEField::ComputeEField(ComputeID c, PatchID pid)
  : ComputePatch(c,pid)
{

	reduction = ReductionMgr::Object()->willSubmit(REDUCTIONS_BASIC);

}
/*			END OF FUNCTION ComputeEField		*/


ComputeEField::~ComputeEField()

{
	delete reduction;
}
/*			END OF FUNCTION ~ComputeEField		*/


void ComputeEField::doForce(CompAtom* p, Results* r) {

  SimParameters *simParams = Node::Object()->simParameters;
  Vector eField = simParams->eField;

  Force *forces = r->f[Results::normal];
  BigReal energy = 0;

  //  Loop through and check each atom
  for (int i=0; i<numAtoms; i++) {
    BigReal charge = p[i].charge;
    forces[i] += charge * eField;
    energy -= charge * ( eField * p[i].position );
  }

  reduction->item(REDUCTION_MISC_ENERGY) += energy;
  reduction->submit();

}
/*			END OF FUNCTION force				*/
