/***************************************************************************/
/*      (C) Copyright 1996,1997 The Board of Trustees of the               */
/*                          University of Illinois                         */
/*                           All Rights Reserved                           */
/***************************************************************************/
/***************************************************************************
 * DESCRIPTION:  Forwards atoms to master node for force evaluation.
 *
 ***************************************************************************/

#ifndef COMPUTEFREEENERGY_H
#define COMPUTEFREEENERGY_H

#include "ComputeGlobalMaster.h"
#include "NamdTypes.h"

class ComputeGlobalConfigMsg;
class ComputeGlobalDataMsg;
class ComputeGlobalResultsMsg;
class ComputeGlobalMaster;
class ComputeMgr;

class ComputeFreeEnergy : public ComputeGlobalMaster {
private:
  friend class ComputeGlobal;
  ComputeFreeEnergy(ComputeGlobal *);
  ~ComputeFreeEnergy();
  virtual void initialize();
  virtual void calculate();
};

#endif
/***************************************************************************
 * RCS INFORMATION:
 *
 *	$RCSfile: ComputeFreeEnergy.h,v $
 *	$Author: jim $	$Locker:  $		$State: Exp $
 *	$Revision: 1.1 $	$Date: 1998/02/10 06:45:09 $
 *
 ***************************************************************************
 * REVISION HISTORY:
 *
 * $Log: ComputeFreeEnergy.h,v $
 * Revision 1.1  1998/02/10 06:45:09  jim
 * Added class ComputeFreeEnergy.
 *
 *
 *
 *
 ***************************************************************************/

