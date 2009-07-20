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
  CkGroupID workDistrib;
  CkGroupID patchMgr;
  CkGroupID proxyMgr;
  CkGroupID computeMgr;
  CkGroupID computePmeMgr;
  CkGroupID computeGridForceNodeMgr;
  CkGroupID computeGridForceMgr;
  //  CkGroupID delegateMgr;
  CkGroupID computeExtMgr;
  CkGroupID reductionMgr;
  CkGroupID collectionMgr;
  CkGroupID broadcastMgr;
  CkGroupID ldbCoordinator;
  CkGroupID sync;
  CkGroupID node;
  
  #if defined(NODEAWARE_PROXY_SPANNINGTREE) && defined(USE_NODEPATCHMGR)
  CkGroupID nodeProxyMgr;
  #endif

};

class BOCclass : public Group {
};

#endif /* BOCGROUP_H */


