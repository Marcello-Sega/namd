#ifndef REDUCTIONMGR_H
#define REDUCTIONMGR_H

#include "NamdTypes.h"

// typedef Vector[3] Tensor;

typedef enum
{
  SUBMIT_KINETIC_ENERGY,
  SUBMIT_BOND_ENERGY,
  SUBMIT_ANGLE_ENERGY,
  SUBMIT_DIHEDRAL_ENERGY,
  SUBMIT_IMPROPER_ENERGY,
  SUBMIT_ELECT_ENERGY,
  SUBMIT_LJ_ENERGY,
  SUBMIT_LONG_RANGE_ENERGY,
  SUBMIT_MAX_RESERVED
} SubmitTag;

typedef enum
{
  REQUIRE_MAX_RESERVED
} RequireTag;

class ReductionMgr
{
public:

  // Vector and Tensor operations can eventually be defined
  // in terms of scalar operations (I hope).
  // For example, reserve the first two tags after a vector
  // and the first eight after a tensor to allow storage.  -JCP

  // (un)register to submit data for reduction
  // raises an error if reductions or broadcasts are active
  register(SubmitTag tag);
  unregister(SubmitTag tag);

  // submit data for reduction
  // more == 1 signals immediate submission of other data
  submit(int seq, SubmitTag tag, BigReal data, int more=0);
  // submit(int seq, SubmitTag tag, Vector data, int more=0);
  // submit(int seq, SubmitTag tag, Tensor data, int more=0);

  // pass on submitting data
  submit(int seq, SubmitTag tag, int more=0);

  // (un)subscribe to broadcast data
  // raises an error if reductions or broadcasts are active
  subscribe(RequireTag tag);
  unsubscribe(RequireTag tag);

  // suspend until this data is ready
  // should be called only from Sequencer thread
  require(int seq, RequireTag tag, BigReal &data);
  // require(int seq, RequireTag tag, Vector &data);
  // require(int seq, RequireTag tag, Tensor &data);

  // pass on requiring data
  require(int seq, RequireTag tag);

  // methods for use by global sequencer

  // raises an error if reductions or broadcasts are active
  request(SubmitTag tag);
  unrequest(SubmitTag tag);

  reduce(int seq, SubmitTag tag, BigReal &data);
  // reduce(int seq, SubmitTag tag, Vector &data);
  // reduce(int seq, SubmitTag tag, Tensor &data);

  reduce(int seq, SubmitTag tag);

  broadcast(int seq, RequireTag tag, BigReal data, int more=0);
  // broadcast(int seq, RequireTag tag, Vector data, int more=0);
  // broadcast(int seq, RequireTag tag, Tensor data, int more=0);

private:

};

#endif
