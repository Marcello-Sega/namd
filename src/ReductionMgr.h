//-*-c++-*-
#ifndef REDUCTIONMGR_H
#define REDUCTIONMGR_H

#include "charm++.h"

#include "main.h"
#include "NamdTypes.h"
#include "BOCgroup.h"
#include "ProcessorPrivate.h"

typedef enum
{
  REDUCTION_ANGLE_ENERGY,
  REDUCTION_BOND_ENERGY,
  REDUCTION_DIHEDRAL_ENERGY,
  REDUCTION_ELECT_ENERGY,
  REDUCTION_IMPROPER_ENERGY,
  REDUCTION_KINETIC_ENERGY,
  REDUCTION_LJ_ENERGY,
  REDUCTION_BC_ENERGY,
  REDUCTION_VIRIAL,
  REDUCTION_ALT_VIRIAL,
  REDUCTION_SMD_ENERGY,
  REDUCTION_MOMENTUM_X,
  REDUCTION_MOMENTUM_Y,
  REDUCTION_MOMENTUM_Z,
  REDUCTION_ANGULAR_MOMENTUM_X,
  REDUCTION_ANGULAR_MOMENTUM_Y,
  REDUCTION_ANGULAR_MOMENTUM_Z,
  REDUCTION_MAX_RESERVED
} ReductionTag;

struct ReductionMgrData
{
  int sequenceNum;
  int numData[REDUCTION_MAX_RESERVED];	// number of data to expect
  int numEvents;	// number of events to expect
  BigReal tagData[REDUCTION_MAX_RESERVED];	// values in tags
  ReductionMgrData *next;	// a queue! ugly but effective.
  int eventCounter;	// counts number of uses

  // for "waiting"
  CthThread threadNum[REDUCTION_MAX_RESERVED];
  int suspendFlag[REDUCTION_MAX_RESERVED];
};

class ReductionDataMsg : public comm_object {
public:
  int seq;
  BigReal data[REDUCTION_MAX_RESERVED];
};

#define MAX_CHILDREN 2

// ***************** for object
class ReductionMgr : public BOCclass
{
private:
  int nextSequence;
  int numSubscribed[REDUCTION_MAX_RESERVED];
  ReductionMgrData *data;	// sequence queue
  int maxData[REDUCTION_MAX_RESERVED];	// number of data to expect
  int maxEvents;	// number of events to expect

  int myParent;
  int numChildren;
  int myChildren[MAX_CHILDREN];

  int isRoot(void) { return (myParent==(-1))?1:0; }
  int isLeaf(void) { return (numChildren==0)?1:0; }

  ReductionMgrData *createdata();		// make new data
  void remove(int seq);				// delete (remove) a sequence
  ReductionMgrData *find(int seq);		// find the data
  void gotAllData(ReductionMgrData *current);	// done collecting data
  void displayData(ReductionMgrData *current);	// display collected data
  void displayData(ReductionMgrData *current, ReductionTag tag);
						// display collected data

public:

  // Singleton Access method
  inline static ReductionMgr *Object(void) {
    return CpvAccess(ReductionMgr_instance);
  }

  ReductionMgr(InitMsg *);
  ~ReductionMgr();

  // (un)register to submit data for reduction
  // may cause an error if reductions are active
  // ASSUMPTION: nobody should register after data has been collected.
  // (register is a key word)
  void Register(ReductionTag tag);
  void unRegister(ReductionTag tag);

  // submit data for reduction
  void submit(int seq, ReductionTag tag, BigReal data);
  void submit(int seq, ReductionTag tag);

  // methods for use by global sequencer

  // raises an error if reductions or broadcasts are active
  void subscribe(ReductionTag tag);
  void unsubscribe(ReductionTag tag);

  // suspend until this data is ready
  // should be called only from Sequencer thread
  void require(int seq, ReductionTag tag, BigReal &data);

  // for receiving data from other ReductionMgr objects
  void recvReductionData(ReductionDataMsg *msg);
};

#endif


/***************************************************************************
 * RCS INFORMATION:
 *
 *	$RCSfile $
 *	$Author $	$Locker:  $		$State: Exp $
 *	$Revision: 1.1010 $	$Date: 1998/04/06 16:34:08 $
 *
 ***************************************************************************
 * REVISION HISTORY:
 *
 * $Log: ReductionMgr.h,v $
 * Revision 1.1010  1998/04/06 16:34:08  jim
 * Added DPME (single processor only), test mode, and momenta printing.
 *
 * Revision 1.1009  1998/03/03 23:05:27  brunner
 * Changed include files for new simplified Charm++ include file structure.
 *
 * Revision 1.1008  1998/01/05 20:27:16  sergei
 * added REDUCTION_SMD_ENERGY to ReductionTag enum.
 *
 * Revision 1.1007  1997/11/07 20:17:49  milind
 * Made NAMD to run on shared memory machines.
 *
 * Revision 1.1006  1997/09/28 10:19:09  milind
 * Fixed priorities, ReductionMgr etc.
 *
 * Revision 1.1005  1997/03/27 03:16:54  jim
 * Added code to check virial calculation, fixed problems with DPMTA and PBC's.
 *
 * Revision 1.1004  1997/03/19 11:54:54  ari
 * Add Broadcast mechanism.
 * Fixed RCS Log entries on files that did not have Log entries.
 * Added some register variables to Molecule and ComputeNonbondedExcl.C
 *
 ***************************************************************************/
