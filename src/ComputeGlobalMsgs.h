/**
***  Copyright (c) 1995, 1996, 1997, 1998, 1999, 2000 by
***  The Board of Trustees of the University of Illinois.
***  All rights reserved.
**/

#ifndef COMPUTEGLOBALMSGS_H
#define COMPUTEGLOBALMSGS_H

#include "charm++.h"

#include "NamdTypes.h"
#include "ComputeMgr.decl.h"

class ComputeGlobalConfigMsg : public CMessage_ComputeGlobalConfigMsg {
public:
  // data members
  AtomIDList aid;
  AtomIDList gdef;  // group definitions

  // constructor and destructor
  ComputeGlobalConfigMsg(void);
  ~ComputeGlobalConfigMsg(void);

  // pack and unpack functions
  static void* pack(ComputeGlobalConfigMsg *msg);
  static ComputeGlobalConfigMsg* unpack(void *ptr);
};


class ComputeGlobalDataMsg : public CMessage_ComputeGlobalDataMsg {
public:
  // data members
  AtomIDList aid;
  PositionList p;
  PositionList gcom;  // group center of mass

  // constructor and destructor
  ComputeGlobalDataMsg(void);
  ~ComputeGlobalDataMsg(void);

  // pack and unpack functions
  static void* pack(ComputeGlobalDataMsg *msg);
  static ComputeGlobalDataMsg* unpack(void *ptr);
};


class ComputeGlobalResultsMsg : public CMessage_ComputeGlobalResultsMsg {
public:
  // data members
  AtomIDList aid;
  ForceList f;
  ForceList gforce;  // force on group
  int reconfig;
  AtomIDList newaid;
  AtomIDList newgdef;

  // constructor and destructor
  ComputeGlobalResultsMsg(void);
  ~ComputeGlobalResultsMsg(void);

  // pack and unpack functions
  static void* pack(ComputeGlobalResultsMsg *msg);
  static ComputeGlobalResultsMsg* unpack(void *ptr);
};


#endif

