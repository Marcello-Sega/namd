#ifndef COMPUTECONSFORCE_H
#define COMPUTECONSFORCE_H

#include "ComputePatch.h"

class ComputeConsForce : public ComputePatch
{
public:
  ComputeConsForce(ComputeID, PatchID);
  virtual void doForce(CompAtom*, Results*);
};

#endif
