//-*-c++-*-
/***************************************************************************/
/*         (C) Copyright 1996,1997 The Board of Trustees of the            */
/*                          University of Illinois                         */
/*                           All Rights Reserved                           */
/***************************************************************************/
/***************************************************************************
 * DESCRIPTION:  Messages needed for ComputeDPME operation.
 *
 ***************************************************************************/

#ifndef COMPUTEDPMEMSGS_H
#define COMPUTEDPMEMSGS_H

#include "charm++.h"

#include "NamdTypes.h"
#include "ComputeMgr.decl.h"

#ifdef DPME
#include "dpme2.h"
#else
#define Pme2Particle char
#define PmeVector char
#endif

class ComputeDPMEDataMsg : public CMessage_ComputeDPMEDataMsg {
public:
  // data members
  int node;
  int numParticles;
  Pme2Particle *particles;

  // constructor and destructor
  ComputeDPMEDataMsg(void);
  ~ComputeDPMEDataMsg(void);

  // pack and unpack functions
  static void* pack (ComputeDPMEDataMsg *msg);
  static ComputeDPMEDataMsg* unpack(void *ptr);
};


class ComputeDPMEResultsMsg : public CMessage_ComputeDPMEResultsMsg {
public:
  // data members
  int node;
  int numParticles;
  PmeVector *forces;

  // constructor and destructor
  ComputeDPMEResultsMsg(void);
  ~ComputeDPMEResultsMsg(void);

  // pack and unpack functions
  static void* pack(ComputeDPMEResultsMsg *msg);
  static ComputeDPMEResultsMsg* unpack(void *ptr);
};


#endif

