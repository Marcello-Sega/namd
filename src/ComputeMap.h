//-*-c++-*-
/***************************************************************************/
/*                                                                         */
/*              (C) Copyright 1996 The Board of Trustees of the            */
/*                          University of Illinois                         */
/*                           All Rights Reserved                           */
/*									   */
/***************************************************************************/

/***************************************************************************
 * DESCRIPTION:
 *
 ***************************************************************************/

#ifndef COMPUTEMAP_H
#define COMPUTEMAP_H

#include "NamdTypes.h"
#include "ProcessorPrivate.h"

class Compute;
class ComputeMgr;

enum ComputeType
{
  computeNonbondedSelfType,
  computeNonbondedPairType,
  computeNonbondedExclType,
  computeBondsType,
  computeAnglesType,
  computeDihedralsType,
  computeImpropersType,
#ifdef DPMTA
  computeDPMTAType,
#endif
#ifdef DPME
  computeDPMEType,
#endif
  computeFullDirectType,
  computeGlobalType,
  computeSphericalBCType,
  computeCylindricalBCType,
  computeRestraintsType,
  computeSMDType,
  computeErrorType
};

class ComputeMap
{
public:
  static ComputeMap *Instance();
  inline static ComputeMap *Object() { return CpvAccess(ComputeMap_instance); }

  void checkMap();

  ~ComputeMap(void);

  void registerCompute(ComputeID cid, Compute *c) {
    computeData[cid].compute = c;
    computeData[cid].moveToNode = -1;
  }

  // numComputes() returns the number of compute objects known
  // by the map.
  int numComputes(void);

  // numPatchBased() returns the number of compute objects
  // that are patch-based
  int numPatchBased(void);

  // numAtomBased() returns the number of compute objects
  // that are atom-based
  int numAtomBased(void);

  // isPatchBased(cid) returns true if the compute object
  // is patch based.
  int isPatchBased(ComputeID cid);

  // isAtomBased(cid) returns true if the compute object
  // is atom based.
  int isAtomBased(ComputeID cid);

  // node(cid) returns the node where the compute object currently exists.
  int node(ComputeID cid);

  void setNode(ComputeID cid, NodeID node);

  // newNode(cid,node) sets up map to tell WorkDistrib to send 
  // compute to new node
  NodeID newNode(ComputeID cid);

  void setNewNode(ComputeID cid, NodeID node);

  // numPids(cid) returns the number of patch ids which are registered
  // with this compute object.
  int numPids(ComputeID cid);
  
  // pid(cid,i) returns the i-th patch id registered
  // with the patch.  
  int pid(ComputeID cid, int i);

  // type(cid) returns the compute type of the given ComputeID
  ComputeType type(ComputeID cid);
  int partition(ComputeID cid);
  int numPartitions(ComputeID cid);

  // allocate_cids(n) tells the ComputeMap to set aside
  // room for n compute objects.  I need not use them all.
  int allocateCids(int n);

  // storeCompute(cid,node,maxPids) tells the ComputeMap to store
  // information about the indicated patch, and allocate space
  // for up to maxPids dependents
  ComputeID storeCompute(int node,int maxPids,ComputeType type,
			 int partition=0, int numPartitions=1);

  // newPid(cid,pid) stores the n patch ids associated with
  // compute id cid.
  int newPid(ComputeID cid, int pid, int trans = 13);

  void printComputeMap(void);

  Compute *compute(ComputeID cid) { return (computeData[cid].compute); };

  friend class ComputeMgr;

  struct PatchRec
  {
    PatchID pid;
    int trans;

    PatchRec() : pid(-1), trans(-1) { ; }
  };

  struct ComputeData
  {
    ComputeData() { 
      node = -1; moveToNode = -1; 
      patchBased = false; numPids = 0; numPidsAllocated = 0; 
      pids = NULL; compute = NULL; 
    }
    Compute *compute;
    int node;
    int moveToNode;
    ComputeType type;
    int partition;
    int numPartitions;
    Boolean patchBased;
    int numPids;
    int numPidsAllocated;
    PatchRec *pids;
  };
protected:
  friend class MapDistribMsg;
  friend class ComputeMapDistribMsg;
  void * pack (int *length);
  void unpack (void *in);

  ComputeMap(void);

private:
  int nPatchBased;
  int nAtomBased;
  int nComputes;
  int nAllocated;
  ComputeData *computeData;
};

#endif /* COMPUTEMAP_H */


/***************************************************************************
 * RCS INFORMATION:
 *
 *	$RCSfile: ComputeMap.h,v $
 *	$Author: brunner $	$Locker:  $		$State: Exp $
 *	$Revision: 1.1015 $	$Date: 1998/07/03 20:09:52 $
 *
 ***************************************************************************
 * REVISION HISTORY:
 *
 * $Log: ComputeMap.h,v $
 * Revision 1.1015  1998/07/03 20:09:52  brunner
 * Self-compute spliting creation changes.  I hope this works.
 *
 * Revision 1.1014  1998/04/06 16:34:04  jim
 * Added DPME (single processor only), test mode, and momenta printing.
 *
 * Revision 1.1013  1998/01/15 04:58:47  jim
 * Corrected "friend foo" to "friend class foo".
 *
 * Revision 1.1012  1998/01/05 20:23:39  sergei
 * added  computeSMDType to enum ComputeType for SMD
 *
 * Revision 1.1011  1997/12/19 23:48:48  jim
 * Added Tcl interface for calculating forces.
 *
 * Revision 1.1010  1997/11/07 20:17:37  milind
 * Made NAMD to run on shared memory machines.
 *
 * Revision 1.1009  1997/07/09 21:26:40  milind
 * Ported NAMD2 to SP3. The SP specific code is within #ifdef SP2
 * and #endif's.
 *
 * Revision 1.1008  1997/04/22 04:25:56  jim
 * Added atomic restraints (harmonic constraints) via ComputeRestraints class.
 *
 * Revision 1.1007  1997/04/10 09:13:50  ari
 * Final debugging for compute migration / proxy creation for load balancing.
 * Lots of debug code added, mostly turned off now.
 * Fixed bug in PositionBox when Patch had no dependencies.
 * Eliminated use of cout and misuse of iout in numerous places.
 *                                            Ari & Jim
 *
 * Revision 1.1006  1997/04/08 07:08:14  ari
 * Modification for dynamic loadbalancing - moving computes
 * Still bug in new computes or usage of proxies/homepatches.
 * Works if ldbStrategy is none as before.
 *
 * Revision 1.1005  1997/03/27 20:25:41  brunner
 * Changes for LdbCoordinator, the load balance control BOC
 *
 * Revision 1.1004  1997/03/20 23:53:39  ari
 * Some changes for comments. Copyright date additions.
 * Hooks for base level update of Compute objects from ComputeMap
 * by ComputeMgr.  Useful for new compute migration functionality.
 *
 * Revision 1.1003  1997/03/19 05:49:54  jim
 * Added ComputeSphericalBC, cleaned up make dependencies.
 *
 * Revision 1.1002  1997/02/11 18:51:42  ari
 * Modified with #ifdef DPMTA to safely eliminate DPMTA codes
 * fixed non-buffering of migration msgs
 * Migration works on multiple processors
 *
 * Revision 1.1001  1997/02/07 22:52:15  jim
 * Eliminated use of nAtomBased and uninitialized memory reads.
 *
 * Revision 1.1000  1997/02/06 15:58:03  ari
 * Resetting CVS to merge branches back into the main trunk.
 * We will stick to main trunk development as suggested by CVS manual.
 * We will set up tags to track fixed points of development/release
 * as suggested by CVS manual - all praise the CVS manual.
 *
 * Revision 1.779  1997/02/06 15:53:02  ari
 * Updating Revision Line, getting rid of branches
 *
 * Revision 1.778.2.1  1997/02/06 02:35:18  jim
 * Implemented periodic boundary conditions - may not work with
 * atom migration yet, but doesn't seem to alter calculation,
 * appears to work correctly when turned on.
 * NamdState chdir's to same directory as config file in argument.
 *
 * Revision 1.778  1997/01/28 00:30:15  ari
 * internal release uplevel to 1.778
 *
 * Revision 1.777.2.1  1997/01/27 22:45:02  ari
 * Basic Atom Migration Code added.
 * Added correct magic first line to .h files for xemacs to go to C++ mode.
 * Compiles and runs without migration turned on.
 *
 * Revision 1.777  1997/01/17 19:35:49  ari
 * Internal CVS leveling release.  Start development code work
 * at 1.777.1.1.
 *
 * Revision 1.7  1996/12/12 08:57:17  jim
 * added MapDistribMsg packing / unpacking routines
 *
 * Revision 1.6  1996/11/30 00:32:52  jim
 * added ComputeMgr friend and storage of ComputeType
 *
 * Revision 1.5  1996/11/01 21:20:45  ari
 * *** empty log message ***
 *
 * Revision 1.4  1996/10/29 23:35:27  ari
 * *** empty log message ***
 *
 * Revision 1.3  1996/10/16 08:22:39  ari
 * *** empty log message ***
 *
 * Revision 1.2  1996/08/16 21:41:11  brunner
 * *** empty log message ***
 *
 * Revision 1.1  1996/08/16 20:43:53  brunner
 * Initial revision
 *
 * Revision 1.2  1996/08/03 20:08:09  brunner
 * *** empty log message ***
 *
 * Revision 1.1  1996/07/16 20:47:43  brunner
 * Initial revision
 *
 *
 ***************************************************************************/
