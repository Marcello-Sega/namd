//-*-c++-*-
#ifndef REDUCTIONMGR_H
#define REDUCTIONMGR_H

#include "ckdefs.h"
#include "chare.h"
#include "c++interface.h"
#include "main.h"
#include "NamdTypes.h"
#include "BOCgroup.h"

// debug code to determine if I should panic
// 1 means PANIC checking is enabled.  0 means no checking.
#ifndef PANIC
#define PANIC	0
#endif

// typedef Vector[3] Tensor;

typedef enum
{
  REDUCTION_ANGLE_ENERGY,
  REDUCTION_BOND_ENERGY,
  REDUCTION_DIHEDRAL_ENERGY,
  REDUCTION_ELECT_ENERGY,
  REDUCTION_IMPROPER_ENERGY,
  REDUCTION_KINETIC_ENERGY,
  REDUCTION_LJ_ENERGY,
  // REDUCTION_LONG_RANGE_ENERGY,
  REDUCTION_MAX_RESERVED
} ReductionTag;

struct ReductionMgrData
{
  int sequenceNum;
  #if PANIC > 0
  int dataToSend;			// number of tags to send (not full)
  #endif
  int numData[REDUCTION_MAX_RESERVED];	// number of data to expect
  int numEvents;	// number of events to expect
  BigReal tagData[REDUCTION_MAX_RESERVED];	// values in tags
  ReductionMgrData *next;	// a queue! ugly but effective.
  int eventCounter;	// counts number of uses

  // for "waiting"
  CthThread threadNum[REDUCTION_MAX_RESERVED];
  int suspendFlag[REDUCTION_MAX_RESERVED];
};

// ***************** for Charm messages
class RegisterReductionMsg : public comm_object
{
  public:
    NodeID node;
    PatchID patch;
};

class UnregisterReductionMsg : public comm_object
{
  public:
    NodeID node;
    PatchID patch;
};

class ReductionDataMsg : public comm_object {
public:
  int seq;
  ReductionTag tag;
  BigReal data;
};


// ***************** for object
class ReductionMgr : public BOCclass
{
private:
  static ReductionMgr *_instance;

  #if PANIC > 0
  int panicMode;
  char tagString[REDUCTION_MAX_RESERVED+1][40];
  #endif

  int nextSequence;
  int numSubscribed[REDUCTION_MAX_RESERVED];
  ReductionMgrData *data;	// sequence queue
  int maxData[REDUCTION_MAX_RESERVED];	// number of data to expect
  int maxEvents;	// number of events to expect

  ReductionMgrData *createdata();		// make new data
  void remove(int seq);				// delete (remove) a sequence
  ReductionMgrData *find(int seq);		// find the data
  void gotAllData(ReductionMgrData *current);	// done collecting data
  void displayData(ReductionMgrData *current);	// display collected data
  void displayData(ReductionMgrData *current, ReductionTag tag);
						// display collected data

public:

  // Singleton Access method
  inline static ReductionMgr *Object() {return _instance;}

  // Vector and Tensor operations can eventually be defined
  // in terms of scalar operations (I hope).
  // For example, reserve the first two tags after a vector
  // and the first eight after a tensor to allow storage.  -JCP

  ReductionMgr(InitMsg *);
  ~ReductionMgr();

  // (un)register to submit data for reduction
  // may cause an error if reductions are active
  // ASSUMPTION: nobody should register after data has been collected.
  // (register is a key word)
  void Register(ReductionTag tag);
  void unRegister(ReductionTag tag);

  // submit data for reduction
  // more == 1 signals immediate submission of other data
  void submit(int seq, ReductionTag tag, BigReal data, int more=0);
  // void submit(int seq, ReductionTag tag, Vector data, int more=0);
  // void submit(int seq, ReductionTag tag, Tensor data, int more=0);

  // pass on submitting data
  void submit(int seq, ReductionTag tag, int more=0);

  // methods for use by global sequencer

  // raises an error if reductions or broadcasts are active
  void subscribe(ReductionTag tag);
  void unsubscribe(ReductionTag tag);

  // suspend until this data is ready
  // should be called only from Sequencer thread
  void require(int seq, ReductionTag tag, BigReal &data);
  // void require(int seq, ReductionTag tag, Vector &data);
  // void require(int seq, ReductionTag tag, Tensor &data);
  // void require(int seq, ReductionTag tag); // pass on requiring data

  // for receiving data from other ReductionMgr objects
  void recvReductionData(ReductionDataMsg *msg);
};

#endif
