/***************************************************************************/
/*      (C) Copyright 1996,1997 The Board of Trustees of the               */
/*                          University of Illinois                         */
/*                           All Rights Reserved                           */
/***************************************************************************/
/***************************************************************************
 * DESCRIPTION:  Forwards atoms to master node for force evaluation.
 *
 ***************************************************************************/

#ifndef COMPUTEGLOBALMISC_H
#define COMPUTEGLOBALMISC_H

#include "ComputeHomePatches.h"
#include "ComputeGlobalEasy.h"
#include "NamdTypes.h"

class ComputeGlobalConfigMsg;
class ComputeGlobalDataMsg;
class ComputeGlobalResultsMsg;
class ComputeGlobalMaster;
class ComputeMgr;

class ComputeMisc : ComputeGlobalEasy {
protected:
  friend class ComputeGlobal;
  ComputeMisc(ComputeGlobal *);
  virtual ~ComputeMisc();

  virtual void easy_init(const char *);
  virtual void easy_calc(void);

};

#endif
/***************************************************************************
 * RCS INFORMATION:
 *
 *	$RCSfile: ComputeMisc.h,v $
 *	$Author: jim $	$Locker:  $		$State: Exp $
 *	$Revision: 1.1 $	$Date: 1999/06/03 16:50:08 $
 *
 ***************************************************************************
 * REVISION HISTORY:
 *
 * $Log: ComputeMisc.h,v $
 * Revision 1.1  1999/06/03 16:50:08  jim
 * Added simplified interface to ComputeGlobal mechanism.
 *
 *
 *
 ***************************************************************************/

