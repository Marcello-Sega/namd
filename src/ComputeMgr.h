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

#ifndef COMPUTEMGR_H
#define COMPUTEMGR_H

#include "charm++.h"
#include "main.h"
#include "new.h"

#include "NamdTypes.h"
#include "BOCgroup.h"

#include "ResizeArray.h"

class Compute;
class ComputeMap;
class CkQdMsg;

class ComputeGlobal;
class ComputeGlobalConfigMsg;
class ComputeGlobalDataMsg;
class ComputeGlobalResultsMsg;

class ComputeDPME;
class ComputeDPMEDataMsg;
class ComputeDPMEResultsMsg;

class ComputePme;
class ComputePmeDataMsg;
class ComputePmeResultsMsg;

class ComputeMgr : public BOCclass
{
public:

  ComputeMgr();
  ~ComputeMgr();
  void createComputes(ComputeMap *map);
  void updateComputes(int,int);
  void updateComputes2(CkQdMsg *);
  void updateComputes3();
  void updateLocalComputes();
  void updateLocalComputes2(CkQdMsg *);
  void updateLocalComputes3();
  void doneUpdateLocalComputes();

  void sendComputeGlobalConfig(ComputeGlobalConfigMsg *);
  void recvComputeGlobalConfig(ComputeGlobalConfigMsg *);
  void sendComputeGlobalData(ComputeGlobalDataMsg *);
  void recvComputeGlobalData(ComputeGlobalDataMsg *);
  void sendComputeGlobalResults(ComputeGlobalResultsMsg *);
  void recvComputeGlobalResults(ComputeGlobalResultsMsg *);

  void sendComputeDPMEData(ComputeDPMEDataMsg *);
  void recvComputeDPMEData(ComputeDPMEDataMsg *);
  void sendComputeDPMEResults(ComputeDPMEResultsMsg *, int);
  void recvComputeDPMEResults(ComputeDPMEResultsMsg *);

  void sendComputePmeData(ComputePmeDataMsg *);
  void recvComputePmeData(ComputePmeDataMsg *);
  void sendComputePmeResults(ComputePmeResultsMsg *, int);
  void recvComputePmeResults(ComputePmeResultsMsg *);

private:
  void createCompute(ComputeID, ComputeMap *);
  int numNonbondedSelf;
  int numNonbondedPair;

  ComputeGlobal *computeGlobalObject;
  ComputeDPME *computeDPMEObject;
  ComputePme *computePmeObject;

  int updateComputesCount;
  int updateComputesReturnEP;
  int updateComputesReturnChareID;

  int *computeFlag;

  class ComputeElem {
  public:
    ComputeID   cid;
    Compute *c;

    int operator<(ComputeElem e) { return (cid < e.cid); }
    int operator==(ComputeElem e) { return (cid == e.cid); }

    ComputeElem(ComputeID id=-1, Compute *compute=NULL) : 
      cid(id), c(compute) {};
    ~ComputeElem() { };
    ComputeElem& operator=(const ComputeElem &e) { 
      cid = e.cid; c = e.c;  // Do not delete c!  This op only used to shuffle
                             // we delete the c here only when the Compute is 
		             // moved off!
      return(*this);
    };
  };

  int numComputes;

  typedef ResizeArray<ComputeElem> ComputeList;
  typedef ResizeArray<int> ComputeIndex;

  // global patch number to local patch table conversion table
  ComputeIndex computeIndex;

  // an array of compute pointers residing on this node
  ComputeList computeList;

  int workDistribGroup;
};

#endif /* COMPUTEMGR_H */

