#ifndef REDUCTIONMGR_H
#define REDUCTIONMGR_H

#include "NamdTypes.h"

// typedef Vector[3] Tensor;

typedef enum
{
  REDUCTION_KINETIC_ENERGY,
  REDUCTION_BOND_ENERGY,
  REDUCTION_ANGLE_ENERGY,
  REDUCTION_DIHEDRAL_ENERGY,
  REDUCTION_IMPROPER_ENERGY,
  REDUCTION_ELECT_ENERGY,
  REDUCTION_LJ_ENERGY,
  REDUCTION_LONG_RANGE_ENERGY,
  REDUCTION_MAX_RESERVED
} ReductionTag;

class ReductionMgr
{
public:

  // Vector and Tensor operations can eventually be defined
  // in terms of scalar operations (I hope).
  // For example, reserve the first two tags after a vector
  // and the first eight after a tensor to allow storage.  -JCP

  // (un)register to submit data for reduction
  // may cause an error if reductions are active
  register(ReductionTag tag);
  unregister(ReductionTag tag);

  // submit data for reduction
  // more == 1 signals immediate submission of other data
  submit(int seq, ReductionTag tag, BigReal data, int more=0);
  submit(int seq, ReductionTag tag, Vector data, int more=0);
  submit(int seq, ReductionTag tag, Tensor data, int more=0);

  // pass on submitting data
  submit(int seq, ReductionTag tag, int more=0);

  // methods for use by global sequencer

  // raises an error if reductions or broadcasts are active
  subscribe(ReductionTag tag);
  unsubscribe(ReductionTag tag);

  // suspend until this data is ready
  // should be called only from Sequencer thread
  require(int seq, ReductionTag tag, BigReal &data);
  require(int seq, ReductionTag tag, Vector &data);
  require(int seq, ReductionTag tag, Tensor &data);
  require(int seq, ReductionTag tag); // pass on requiring data

private:

};

#endif
