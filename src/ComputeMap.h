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

class Compute;
class ComputeMgr;

class ComputeMap
{
public:
  static ComputeMap *Instance();
  inline static ComputeMap *Object() { return _instance; }

  ~ComputeMap(void);

  void registerCompute(ComputeID cid, Compute *c) {
    computeData[cid].compute = c;
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
  int isPatchBased(int cid);

  // isAtomBased(cid) returns true if the compute object
  // is atom based.
  int isAtomBased(int cid);

  // node(cid) returns the node where the compute object currently exists.
  int node(int cid);

  // numPids(cid) returns the number of patch ids which are registered
  // with this compute object.
  int numPids(int cid);
  
  // pid(cid,i) returns the i-th patch id registered
  // with the patch.  
  int pid(int cid, int i);

  // allocate_cids(n) tells the ComputeMap to set aside
  // room for n compute objects.  I need not use them all.
  int allocateCids(int n);

  // storeCompute(cid,node,maxPids) tells the ComputeMap to store
  // information about the indicated patch, and allocate space
  // for up to maxPids dependents
  ComputeID storeCompute(int node,int maxPids,ComputeType type);

  // newPid(cid,pid) stores the n patch ids associated with
  // compute id cid.
  int newPid(int cid, int pid, int trans = 13);

  void printComputeMap(void);

  Compute *compute(ComputeID cid) { return (computeData[cid].compute); };

  friend ComputeMgr;

protected:
  friend MapDistribMsg;
  void * pack (int *length);
  void unpack (void *in);

  ComputeMap(void);

private:
  static ComputeMap *_instance;

  struct PatchRec
  {
    PatchID pid;
    int trans;

    PatchRec() : pid(-1), trans(-1) { ; }
  };

  struct ComputeData
  {
    Compute *compute;
    int node;
    ComputeType type;
    Boolean patchBased;
    int numPids;
    int numPidsAllocated;
    PatchRec *pids;
  };

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
 *	$Author: jim $	$Locker:  $		$State: Exp $
 *	$Revision: 1.1001 $	$Date: 1997/02/07 22:52:15 $
 *
 ***************************************************************************
 * REVISION HISTORY:
 *
 * $Log: ComputeMap.h,v $
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
