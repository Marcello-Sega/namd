/***************************************************************************/
/*      (C) Copyright 1996,1997 The Board of Trustees of the               */
/*                          University of Illinois                         */
/*                           All Rights Reserved                           */
/***************************************************************************/
/***************************************************************************
 * DESCRIPTION:  Computes electrostatics using efficient all atoms^2 calc
 *
 ***************************************************************************/

#ifndef COMPUTEDPME_H
#define COMPUTEDPME_H

#include "ComputeHomePatches.h"
#include "NamdTypes.h"

#ifdef DPME

class ComputeDPMEDataMsg;
class ComputeDPMEResultsMsg;
class ComputeDPMEMaster;
class ComputeMgr;

class ComputeDPME : public ComputeHomePatches {
public:
  ComputeDPME(ComputeID c, ComputeMgr *m);
  virtual ~ComputeDPME();
  void doWork();
  void recvData(ComputeDPMEDataMsg *);
  void recvResults(ComputeDPMEResultsMsg *);

  ComputeMgr *comm;
  int getMasterNode(void) { return masterNode; }

 private:
  ComputeDPMEMaster *master;
  int masterNode;
  int numLocalAtoms;

};

#endif
#endif

