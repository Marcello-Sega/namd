
#ifndef COMPUTEMSM_H
#define COMPUTEMSM_H

#include "Lattice.h"
#include "ComputeMsmMgr.decl.h"
#include "ComputeHomePatches.h"
#include "NamdTypes.h"


class ComputeMsmMgr;
class SubmitReduction;

class MsmInitMsg : public CMessage_MsmInitMsg {
  public:
    ScaledPosition smin, smax;  // needs the extreme positions
};


class ComputeMsm : public ComputeHomePatches {
public:
  ComputeMsm(ComputeID c);
  virtual ~ComputeMsm();
  void doWork();
  void saveResults(/* int n, const Force [], double self_energy */);

  void setMgr(ComputeMsmMgr *mgr) { myMgr = mgr; }

private:
  double qscaling;  // charge scaling constant
  SubmitReduction *reduction;

  ComputeMsmMgr *myMgr;
};

#if 0
struct MsmData {
  int ispx, ispy, ispz;
  void pup(PUP::er &p);  // for parameter marshalling
  void print();          // for debugging
};
#endif


#endif // COMPUTEMSM_H

