//-*-c++-*-
/***************************************************************************/
/*                                                                         */
/*              (C) Copyright 1996 The Board of Trustees of the            */
/*                          University of Illinois                         */
/*                           All Rights Reserved                           */
/*                                                                         */
/***************************************************************************/

/***************************************************************************
 * DESCRIPTION:
 *
 ***************************************************************************/

#ifndef BOCGROUP_H
#define BOCGROUP_H

#include "charm++.h"

class BOCgroup {
public:
  int workDistrib;
  int patchMgr;
  int proxyMgr;
  int computeMgr;
  int reductionMgr;
  int collectionMgr;
  int broadcastMgr;
  int ldbCoordinator;
  int node;
};

class BOCclass : public Group {
};

#endif /* BOCGROUP_H */


