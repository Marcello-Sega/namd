/***************************************************************************/
/*      (C) Copyright 1996,1997 The Board of Trustees of the               */
/*                          University of Illinois                         */
/*                           All Rights Reserved                           */
/***************************************************************************/
/***************************************************************************
 * DESCRIPTION:  Forwards atoms to master node for force evaluation.
 *
 ***************************************************************************/

#ifndef COMPUTEGLOBALMASTER_H
#define COMPUTEGLOBALMASTER_H

#include "ComputeHomePatches.h"
#include "NamdTypes.h"

class ComputeGlobalConfigMsg;
class ComputeGlobalDataMsg;
class ComputeGlobalResultsMsg;
class ComputeGlobalMaster;
class ComputeMgr;

class ComputeGlobalMaster {
protected:
  friend class ComputeGlobal;
  ComputeGlobal *host;
  ComputeGlobalMaster(ComputeGlobal *);
  virtual ~ComputeGlobalMaster();
  void recvData(ComputeGlobalDataMsg *);
  int numWorkingPes;
  int msgcount;
  virtual void initialize();
  int initialized;
  void storedata(ComputeGlobalDataMsg *);
  void cleardata();
  AtomIDList aid;
  PositionList p;
  PositionList gcom;  // group centers of mass
  void storedefs(AtomIDList newgdef);
  AtomIDList gdef;  // group definitions
  ResizeArray<BigReal> gmass;  // group masses
  virtual void calculate();
};

#endif

