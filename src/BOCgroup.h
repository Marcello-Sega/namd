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

#include "ckdefs.h"
#include "chare.h"
#include "c++interface.h"

class BOCgroup {
public:
  int workDistrib;
  int patchMgr;
  int proxyMgr;
  int computeMgr;
  int reductionMgr;
  int collectionMgr;
  int node;
};

class BOCclass : public groupmember {
protected:
    static BOCgroup group;
};

#endif /* BOCGROUP_H */
