//-*-c++-*-
#ifndef REDUCTIONMGR_H
#define REDUCTIONMGR_H

#include "charm++.h"

#include "main.h"
#include "NamdTypes.h"
#include "BOCgroup.h"
#include "ProcessorPrivate.h"
#include "ReductionMgr.decl.h"

typedef enum
{
 // energy
  REDUCTION_ANGLE_ENERGY,
  REDUCTION_BOND_ENERGY,
  REDUCTION_DIHEDRAL_ENERGY,
  REDUCTION_ELECT_ENERGY,
  REDUCTION_ELECT_ENERGY_SLOW,
  REDUCTION_IMPROPER_ENERGY,
  REDUCTION_KINETIC_ENERGY,
  REDUCTION_INT_KINETIC_ENERGY,
  REDUCTION_LJ_ENERGY,
  REDUCTION_BC_ENERGY,
  REDUCTION_SMD_ENERGY,
  REDUCTION_MISC_ENERGY,
 // pressure
  REDUCTION_VIRIAL_NORMAL_X,
  REDUCTION_VIRIAL_NORMAL_Y,
  REDUCTION_VIRIAL_NORMAL_Z,
  REDUCTION_VIRIAL_NBOND_X,
  REDUCTION_VIRIAL_NBOND_Y,
  REDUCTION_VIRIAL_NBOND_Z,
  REDUCTION_VIRIAL_SLOW_X,
  REDUCTION_VIRIAL_SLOW_Y,
  REDUCTION_VIRIAL_SLOW_Z,
  REDUCTION_ALT_VIRIAL_NORMAL_X,
  REDUCTION_ALT_VIRIAL_NORMAL_Y,
  REDUCTION_ALT_VIRIAL_NORMAL_Z,
  REDUCTION_ALT_VIRIAL_NBOND_X,
  REDUCTION_ALT_VIRIAL_NBOND_Y,
  REDUCTION_ALT_VIRIAL_NBOND_Z,
  REDUCTION_ALT_VIRIAL_SLOW_X,
  REDUCTION_ALT_VIRIAL_SLOW_Y,
  REDUCTION_ALT_VIRIAL_SLOW_Z,
  REDUCTION_INT_VIRIAL_NORMAL_X,
  REDUCTION_INT_VIRIAL_NORMAL_Y,
  REDUCTION_INT_VIRIAL_NORMAL_Z,
  REDUCTION_INT_VIRIAL_NBOND_X,
  REDUCTION_INT_VIRIAL_NBOND_Y,
  REDUCTION_INT_VIRIAL_NBOND_Z,
  REDUCTION_INT_VIRIAL_SLOW_X,
  REDUCTION_INT_VIRIAL_SLOW_Y,
  REDUCTION_INT_VIRIAL_SLOW_Z,
 // momentum
  REDUCTION_MOMENTUM_X,
  REDUCTION_MOMENTUM_Y,
  REDUCTION_MOMENTUM_Z,
  REDUCTION_ANGULAR_MOMENTUM_X,
  REDUCTION_ANGULAR_MOMENTUM_Y,
  REDUCTION_ANGULAR_MOMENTUM_Z,
 // checksum
  REDUCTION_ATOM_CHECKSUM,
  REDUCTION_COMPUTE_CHECKSUM,
  REDUCTION_BOND_CHECKSUM,
  REDUCTION_ANGLE_CHECKSUM,
  REDUCTION_DIHEDRAL_CHECKSUM,
  REDUCTION_IMPROPER_CHECKSUM,
  REDUCTION_EXCLUSION_CHECKSUM,
  REDUCTION_MARGIN_VIOLATIONS,
 // semaphore (must be last)
  REDUCTION_MAX_RESERVED
} ReductionTag;

// Later this can be dynamic
enum {
  REDUCTIONS_BASIC,
  REDUCTIONS_USER1,
  REDUCTIONS_USER2,
 // semaphore (must be last)
  REDUCTION_MAX_SET_ID
};

class ReductionRegisterMsg;
class ReductionSubmitMsg;
class ReductionSet;
class SubmitReduction;
class RequireReduction;

// Top level class
class ReductionMgr : public BOCclass
{
private:
  friend class SubmitReduction;
  friend class RequireReduction;

  ReductionSet * (reductionSets[REDUCTION_MAX_SET_ID]);

  int myParent;  // parent node or -1 if none
  int firstChild, lastChild;  // firstChild <= children < lastChild
  int isMyChild(int nodeID) const {
    return ( nodeID >= firstChild && nodeID < lastChild );
  }
  int isRoot(void) const { return ( myParent == -1 ); }

  ReductionSet* getSet(int setID);
  void delSet(int setID);

  void mergeAndDeliver(
	ReductionSet *set, int seqNum, const BigReal *newData, int size);

  void submit(SubmitReduction*);
  void remove(SubmitReduction*);

  void require(RequireReduction*);
  void remove(RequireReduction*);

public:

  // Singleton Access method
  inline static ReductionMgr *Object(void) {
    return CpvAccess(ReductionMgr_instance);
  }

  ReductionMgr();
  ~ReductionMgr();

  // client interface
  SubmitReduction* willSubmit(int setID);
  RequireReduction* willRequire(int setID);

  // message entry points
  void remoteRegister(ReductionRegisterMsg *msg);
  void remoteUnregister(ReductionRegisterMsg *msg);
  void remoteSubmit(ReductionSubmitMsg *msg);

};

// Client handle for submissions
class SubmitReduction {
private:
  friend class ReductionMgr;
  int reductionSetID;
  int sequenceNumber;
  ReductionMgr *master;
  int dataSize;
  BigReal *data;
  SubmitReduction(void) { dataSize = 0; data = 0; }
public:
  BigReal& item(int i) {
    if ( i >= dataSize ) {
      int oldSize = dataSize;
      BigReal *oldData = data;
      dataSize = i+1;
      data = new BigReal[dataSize];
      int j = 0;
      for ( ; j < oldSize; ++j ) data[j] = oldData[j];
      for ( ; j < dataSize; ++j ) data[j] = 0.;
      delete [] oldData;
    }
    return data[i];
  }
  void submit(void) {
    master->submit(this);
    ++sequenceNumber;
    for ( int i = 0; i < dataSize; ++i ) { data[i] = 0; }
  }
  ~SubmitReduction(void) { delete [] data; master->remove(this); }
};

// Client handle for requires
class RequireReduction {
private:
  friend class ReductionMgr;
  int reductionSetID;
  int sequenceNumber;
  ReductionMgr *master;
  int dataSize;
  BigReal *data;
  RequireReduction(void) { dataSize = 0; data = 0; }
public:
  BigReal item(int i) const { return ( i < dataSize ? data[i] : 0 ); }
  void require(void) {
    master->require(this);
    ++sequenceNumber;
  }
  ~RequireReduction(void) { delete [] data; master->remove(this); }
};


#endif

