
#ifndef COMPUTE_IMD_H__
#define COMPUTE_IMD_H__

#include "ComputeGlobalMaster.h"

class ComputeGlobalConfigMsg;
class ComputeGlobalResultsMsg;

class ComputeIMD : public ComputeGlobalMaster {
friend class ComputeGlobal;

private:
  ComputeIMD(ComputeGlobal *);
  ~ComputeIMD();

  virtual void initialize();
  virtual void calculate();

  ComputeGlobalConfigMsg *configMsg;
  ComputeGlobalResultsMsg *resultsMsg;

protected:
  // My socket handle
  void *sock;

  // Simple function for getting MDComm-style forces from VMD
  int get_vmd_forces();

  // These are the data structures for MDCOMM-style sending of indices
  // and forces from VMD to NAMD.
  int num_vmd_atoms;
  int *vmd_atoms;
  float *vmd_forces;
};

#endif

