/***************************************************************************/
/*      (C) Copyright 1996,1997 The Board of Trustees of the               */
/*                          University of Illinois                         */
/*                           All Rights Reserved                           */
/***************************************************************************/
/***************************************************************************
 * DESCRIPTION:  Forwards atoms to master node for force evaluation.
 *
 ***************************************************************************/

#ifndef COMPUTEGLOBAL_H
#define COMPUTEGLOBAL_H

#include "ComputeHomePatches.h"
#include "NamdTypes.h"

class ComputeGlobalConfigMsg;
class ComputeGlobalDataMsg;
class ComputeGlobalResultsMsg;
class ComputeGlobalMaster;
class ComputeMgr;

class ComputeGlobal : public ComputeHomePatches {
public:
  ComputeGlobal(ComputeID, ComputeMgr*);
  virtual ~ComputeGlobal();
  void doWork();
  void recvConfig(ComputeGlobalConfigMsg *);
  void recvData(ComputeGlobalDataMsg *);
  void recvResults(ComputeGlobalResultsMsg *);

private:
  friend class ComputeGlobalMaster;
  ComputeGlobalMaster *master;
  ComputeMgr *comm;

  void sendData();
  int configured;
  AtomIDList aid;
};

#endif
/***************************************************************************
 * RCS INFORMATION:
 *
 *	$RCSfile: ComputeGlobal.h,v $
 *	$Author: jim $	$Locker:  $		$State: Exp $
 *	$Revision: 1.2 $	$Date: 1998/01/15 04:58:46 $
 *
 ***************************************************************************
 * REVISION HISTORY:
 *
 * $Log: ComputeGlobal.h,v $
 * Revision 1.2  1998/01/15 04:58:46  jim
 * Corrected "friend foo" to "friend class foo".
 *
 * Revision 1.1  1997/12/19 23:48:46  jim
 * Added Tcl interface for calculating forces.
 *
 *
 ***************************************************************************/

