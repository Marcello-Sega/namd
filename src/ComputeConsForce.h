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

#ifndef COMPUTECONSTORQUE_H
#define COMPUTECONSTORQUE_H

#include "ComputePatch.h"

class ComputeConsTorque : public ComputePatch
{
public:
  ComputeConsTorque(ComputeID, PatchID);
  virtual void doForce(CompAtom*, Results*);
};

#endif
