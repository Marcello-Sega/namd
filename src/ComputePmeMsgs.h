//-*-c++-*-
/***************************************************************************/
/*         (C) Copyright 1996,1997 The Board of Trustees of the            */
/*                          University of Illinois                         */
/*                           All Rights Reserved                           */
/***************************************************************************/
/***************************************************************************
 * DESCRIPTION:  Messages needed for ComputePme operation.
 *
 ***************************************************************************/

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
  int node;
  int numParticles;
  Vector *forces;

  // constructor and destructor
  ComputePmeResultsMsg(void);
  ~ComputePmeResultsMsg(void);

  // pack and unpack functions
  static void* pack(ComputePmeResultsMsg *msg);
  static ComputePmeResultsMsg* unpack(void *ptr);
};


#endif

