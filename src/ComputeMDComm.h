/***************************************************************************/
/*      (C) Copyright 1996,1997 The Board of Trustees of the               */
/*                          University of Illinois                         */
/*                           All Rights Reserved                           */
/***************************************************************************/
/***************************************************************************
 * DESCRIPTION:  Forwards atoms to master node for force evaluation.
 *
 ***************************************************************************/

#ifndef COMPUTEMDCOMM_H
#define COMPUTEMDCOMM_H

#include "ComputeGlobalMaster.h"
#include "NamdTypes.h"
#include "iostream.h"

class ComputeGlobalConfigMsg;
class ComputeGlobalDataMsg;
class ComputeGlobalResultsMsg;
class ComputeGlobalMaster;
class ComputeMgr;

#ifdef MDCOMM
void mdcomm_transfer_vmdForceData(int, int*, float*);
#endif

class ComputeMDComm : public ComputeGlobalMaster {
private:
  friend class ComputeGlobal;
#ifdef MDCOMM
  friend void mdcomm_transfer_vmdForceData(int, int*, float*);
#endif
  ComputeMDComm(ComputeGlobal *);
  ~ComputeMDComm();
  virtual void initialize();
  virtual void calculate();
  ComputeGlobalConfigMsg *configMsg;
  ComputeGlobalResultsMsg *resultsMsg;
protected:
  static AtomIDList vmd_atoms;
  static ForceList vmd_forces;
};

#endif
/***************************************************************************
 * RCS INFORMATION:
 *
 *	$RCSfile: ComputeMDComm.h,v $
 *	$Author: jim $	$Locker:  $		$State: Exp $
 *	$Revision: 1.1 $	$Date: 1998/04/30 04:53:23 $
 *
 ***************************************************************************
 * REVISION HISTORY:
 *
 * $Log: ComputeMDComm.h,v $
 * Revision 1.1  1998/04/30 04:53:23  jim
 * Added forces from MDComm and other improvements to ComputeGlobal.
 *
 *
 *
 ***************************************************************************/

