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

#ifdef DPME

#include "ComputeHomePatches.h"
#include "NamdTypes.h"

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
/***************************************************************************
 * RCS INFORMATION:
 *
 *	$RCSfile: ComputeDPME.h,v $
 *	$Author: jim $	$Locker:  $		$State: Exp $
 *	$Revision: 1.2 $	$Date: 1998/04/10 04:15:57 $
 *
 ***************************************************************************
 * REVISION HISTORY:
 *
 * $Log: ComputeDPME.h,v $
 * Revision 1.2  1998/04/10 04:15:57  jim
 * Finished incorporating DPME.
 *
 * Revision 1.1  1998/04/07 00:52:38  jim
 * Added DPME interface.
 *
 *
 ***************************************************************************/

