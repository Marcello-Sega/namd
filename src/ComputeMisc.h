/**
***  Copyright (c) 1995, 1996, 1997, 1998, 1999, 2000 by
***  The Board of Trustees of the University of Illinois.
***  All rights reserved.
**/

/*
   Forwards atoms to master node for force evaluation.
*/

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

class ComputeMisc : public ComputeGlobalEasy {
protected:
  friend class ComputeGlobal;
  ComputeMisc(ComputeGlobal *);
  virtual ~ComputeMisc();

  virtual void easy_init(const char *);
  virtual void easy_calc(void);

};

#endif

