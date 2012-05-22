
#ifndef COMPUTEMSM_H
#define COMPUTEMSM_H

#include "Lattice.h"
#include "ComputeMsmMgr.decl.h"
#include "ComputeHomePatches.h"
#include "NamdTypes.h"

class SubmitReduction;
typedef Force MsmForce;

class MsmInitMsg : public CMessage_MsmInitMsg {
  public:
    ScaledPosition smin, smax;  // needs the extreme positions
};


class ComputeMsm : public ComputeHomePatches {
public:
  ComputeMsm(ComputeID c);
  virtual ~ComputeMsm();
  void doWork();
  void saveResults(int n, const MsmForce [], double self_energy);

private:
  double qscaling;  // charge scaling constant
  SubmitReduction *reduction;
};

struct MsmData {
  int ispx, ispy, ispz;
  void pup(PUP::er &p);  // for parameter marshalling
  void print();          // for debugging
};


#endif // COMPUTEMSM_H

