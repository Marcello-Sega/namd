/**
***  Copyright (c) 1995, 1996, 1997, 1998, 1999, 2000 by
***  The Board of Trustees of the University of Illinois.
***  All rights reserved.
**/

#ifndef BOCGROUP_H
#define BOCGROUP_H

#include "charm++.h"

class BOCgroup {
public:
  int workDistrib;
  int patchMgr;
  int proxyMgr;
  int computeMgr;
  int computePmeMgr;
  int reductionMgr;
  int collectionMgr;
  int broadcastMgr;
  int ldbCoordinator;
  int sync;
  int node;
};

class BOCclass : public Group {
};

#endif /* BOCGROUP_H */


