/**
***  Copyright (c) 1995, 1996, 1997, 1998, 1999, 2000 by
***  The Board of Trustees of the University of Illinois.
***  All rights reserved.
**/

#ifndef REDUCTIONMGR_H
#define REDUCTIONMGR_H

#include "charm++.h"

#include "main.h"
#include "NamdTypes.h"
#include "BOCgroup.h"
#include "ProcessorPrivate.h"

#define VECTOR(A) A ## _X, A ## _Y, A ## _Z
#define TENSOR(A) A ## _XX, A ## _XY, A ## _XZ, \
                  A ## _YX, A ## _YY, A ## _YZ, \
                  A ## _ZX, A ## _ZY, A ## _ZZ

#define ADD_VECTOR(R,RL,D,DL) \
  R->item( RL ## _X ) += D[ DL ## _X ]; \
  R->item( RL ## _Y ) += D[ DL ## _Y ]; \
  R->item( RL ## _Z ) += D[ DL ## _Z ]

#define ADD_TENSOR(R,RL,D,DL) \
  R->item( RL ## _XX) += D[ DL ## _XX ]; \
  R->item( RL ## _XY) += D[ DL ## _XY ]; \
  R->item( RL ## _XZ) += D[ DL ## _XZ ]; \
  R->item( RL ## _YX) += D[ DL ## _YX ]; \
  R->item( RL ## _YY) += D[ DL ## _YY ]; \
  R->item( RL ## _YZ) += D[ DL ## _YZ ]; \
  R->item( RL ## _ZX) += D[ DL ## _ZX ]; \
  R->item( RL ## _ZY) += D[ DL ## _ZY ]; \
  R->item( RL ## _ZZ) += D[ DL ## _ZZ ]

#define ADD_TENSOR_OBJECT(R,RL,D) \
  R->item( RL ## _XX) += D.xx; \
  R->item( RL ## _XY) += D.xy; \
  R->item( RL ## _XZ) += D.xz; \
  R->item( RL ## _YX) += D.yx; \
  R->item( RL ## _YY) += D.yy; \
  R->item( RL ## _YZ) += D.yz; \
  R->item( RL ## _ZX) += D.zx; \
  R->item( RL ## _ZY) += D.zy; \
  R->item( RL ## _ZZ) += D.zz

#define GET_VECTOR(O,R,A) \
  O.x = R->item( A ## _X ); \
  O.y = R->item( A ## _Y ); \
  O.z = R->item( A ## _Z )

#define GET_TENSOR(O,R,A) \
  O.xx = R->item( A ## _XX); \
  O.xy = R->item( A ## _XY); \
  O.xz = R->item( A ## _XZ); \
  O.yx = R->item( A ## _YX); \
  O.yy = R->item( A ## _YY); \
  O.yz = R->item( A ## _YZ); \
  O.zx = R->item( A ## _ZX); \
  O.zy = R->item( A ## _ZY); \
  O.zz = R->item( A ## _ZZ)

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
  TENSOR(REDUCTION_VIRIAL_NORMAL),
  TENSOR(REDUCTION_VIRIAL_NBOND),
  TENSOR(REDUCTION_VIRIAL_SLOW),
#ifdef ALTVIRIAL
  TENSOR(REDUCTION_ALT_VIRIAL_NORMAL),
  TENSOR(REDUCTION_ALT_VIRIAL_NBOND),
  TENSOR(REDUCTION_ALT_VIRIAL_SLOW),
#endif
  TENSOR(REDUCTION_INT_VIRIAL_NORMAL),
  TENSOR(REDUCTION_INT_VIRIAL_NBOND),
  TENSOR(REDUCTION_INT_VIRIAL_SLOW),
 // momentum
  VECTOR(REDUCTION_MOMENTUM),
  VECTOR(REDUCTION_ANGULAR_MOMENTUM),
 // used for minimization
  REDUCTION_MIN_F_DOT_F,
  REDUCTION_MIN_F_DOT_V,
  REDUCTION_MIN_V_DOT_V,
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

