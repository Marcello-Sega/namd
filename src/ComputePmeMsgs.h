/**
***  Copyright (c) 1995, 1996, 1997, 1998, 1999, 2000 by
***  The Board of Trustees of the University of Illinois.
***  All rights reserved.
**/

#ifndef COMPUTEPMEMSGS_H
#define COMPUTEPMEMSGS_H

#include "charm++.h"

#include "NamdTypes.h"
#include "ComputeMgr.decl.h"

#include "PmeBase.h"

class ComputePmeDataMsg : public CMessage_ComputePmeDataMsg {
public:
  // data members
  int node;
  int numParticles;
  PmeParticle *particles;

  // constructor and destructor
  ComputePmeDataMsg(void);
  ~ComputePmeDataMsg(void);

  // pack and unpack functions
  static void* pack (ComputePmeDataMsg *msg);
  static ComputePmeDataMsg* unpack(void *ptr);
};


class ComputePmeResultsMsg : public CMessage_ComputePmeResultsMsg {
public:
  // data members
  double energy;
  double virial[6];
  int start;
  int q_len;
  double *q_arr;
  // not really data
  int will_delete_array;

  // constructor and destructor
  ComputePmeResultsMsg(void);
  ~ComputePmeResultsMsg(void);

  // pack and unpack functions
  static void* pack(ComputePmeResultsMsg *msg);
  static ComputePmeResultsMsg* unpack(void *ptr);
};


#endif

