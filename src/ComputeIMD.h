/**
***  Copyright (c) 1995, 1996, 1997, 1998, 1999, 2000 by
***  The Board of Trustees of the University of Illinois.
***  All rights reserved.
**/

#ifndef COMPUTE_IMD_H__
#define COMPUTE_IMD_H__

#include "ComputeGlobalMaster.h"

class ComputeGlobalConfigMsg;
class ComputeGlobalResultsMsg;

class ComputeIMD : public ComputeGlobalMaster {
friend class ComputeGlobal;

private:
  ComputeIMD(ComputeGlobal *);
  ~ComputeIMD();

  virtual void initialize();
  virtual void calculate();

  ComputeGlobalConfigMsg *configMsg;
  ComputeGlobalResultsMsg *resultsMsg;

  // Simple function for getting MDComm-style forces from VMD
  int get_vmd_forces();

  // flag for whether to proceed with simulation when there are no connections
  int IMDwait;

  // My socket handle
  void *sock;
  void *clientsock;
};

#endif

