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

  ComputeMgr *comm;

private:
  ComputeGlobalMaster *master;

  void sendData();
  int configured;
  void configure(AtomIDList newaid, AtomIDList newgdef);
  AtomIDList aid;
  AtomIDList gdef;  // definitions of groups
  ResizeArray<BigReal> gmass;  // masses of groups
};

#endif
/***************************************************************************
 * RCS INFORMATION:
 *
 *	$RCSfile: ComputeGlobal.h,v $
 *	$Author: jim $	$Locker:  $		$State: Exp $
 *	$Revision: 1.4 $	$Date: 1998/02/16 00:23:19 $
 *
 ***************************************************************************
 * REVISION HISTORY:
 *
 * $Log: ComputeGlobal.h,v $
 * Revision 1.4  1998/02/16 00:23:19  jim
 * Added atom group centers of mass to Tcl interface.
 *
 * Revision 1.3  1998/02/10 05:35:02  jim
 * Split ComputeGlobal into different classes and files.
 * Switched globalForces and globalForcesTcl to tclForces and tclForcesScript.
 * Added (soon to be used) freeEnergy and freeEnergyConfig.
 *
 * Revision 1.2  1998/01/15 04:58:46  jim
 * Corrected "friend foo" to "friend class foo".
 *
 * Revision 1.1  1997/12/19 23:48:46  jim
 * Added Tcl interface for calculating forces.
 *
 *
 ***************************************************************************/

