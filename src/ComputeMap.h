/***************************************************************************/
/*                                                                         */
/*              (C) Copyright 1996 The Board of Trustees of the            */
/*                          University of Illinois                         */
/*                           All Rights Reserved                           */
/*									   */
/***************************************************************************/

/***************************************************************************
 * RCS INFORMATION:
 *
 *	$RCSfile: ComputeMap.h,v $
 *	$Author: brunner $	$Locker:  $		$State: Exp $
 *	$Revision: 1.1 $	$Date: 1996/08/16 20:43:53 $
 *
 ***************************************************************************
 * DESCRIPTION:
 *
 ***************************************************************************
 * REVISION HISTORY:
 *
 * $Log: ComputeMap.h,v $
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

#ifndef COMPUTEMAP_H
#define COMPUTEMAP_H

#include "NamdTypes.h"

class ComputeMap
{
private:
  struct ComputeData
  {
    int node;
    int patchBased;
    int numPids;
    int numPidsAllocated;
    PatchID *pids;
  };

  int nPatchBased;
  int nAtomBased;
  int nComputes;
  int nAllocated;
  ComputeData *computeData;
  
public:
  ComputeMap(void);

  ~ComputeMap(void);

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
  int newPid(int cid, int pid);

  void printComputeMap(void);

};

#endif /* COMPUTEMAP_H */






