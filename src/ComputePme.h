/***************************************************************************/
/*      (C) Copyright 1996,1997 The Board of Trustees of the               */
/*                          University of Illinois                         */
/*                           All Rights Reserved                           */
/***************************************************************************/
/***************************************************************************
 * DESCRIPTION:  Computes electrostatics using efficient all atoms^2 calc
 *
 ***************************************************************************/

#ifndef COMPUTEPME_H
#define COMPUTEPME_H

#include "ComputeHomePatches.h"
#include "NamdTypes.h"

#ifdef DPME

class ComputePmeDataMsg;
class ComputePmeResultsMsg;
class ComputePmeMaster;
class ComputeMgr;

class ComputePme : public ComputeHomePatches {
public:
  ComputePme(ComputeID c, ComputeMgr *m);
  virtual ~ComputePme();
  void doWork();
  void recvData(ComputePmeDataMsg *);
  void recvResults(ComputePmeResultsMsg *);

  ComputeMgr *comm;
  int getMasterNode(void) { return masterNode; }

 private:
  ComputePmeMaster *master;
  int masterNode;
  int numLocalAtoms;

};

#endif
#endif
/***************************************************************************
 * RCS INFORMATION:
 *
 *	$RCSfile: ComputePme.h,v $
 *	$Author: jim $	$Locker:  $		$State: Exp $
 *	$Revision: 1.1 $	$Date: 1999/06/08 14:52:07 $
 *
 ***************************************************************************
 * REVISION HISTORY:
 *
 * $Log: ComputePme.h,v $
 * Revision 1.1  1999/06/08 14:52:07  jim
 * Incorporated Justin's faster PME code along side DPME.
 *
 * Revision 1.3  1998/04/15 22:13:50  jim
 * Make depends returns same results regardless of DPME, DPMTA, TCL or MDCOMM.
 *
 * Revision 1.2  1998/04/10 04:15:57  jim
 * Finished incorporating DPME.
 *
 * Revision 1.1  1998/04/07 00:52:38  jim
 * Added DPME interface.
 *
 *
 ***************************************************************************/

