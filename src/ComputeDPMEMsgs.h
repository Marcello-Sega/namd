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

#ifdef DPME
#include "dpme2.h"
#else
class Pme2Particle;
class PmeVector;
#endif

class ComputeDPMEDataMsg : public comm_object {
public:
  // data members
  int node;
  int numParticles;
  Pme2Particle *particles;

  // constructor and destructor
  ComputeDPMEDataMsg(void);
  ~ComputeDPMEDataMsg(void);

  // standard new overload for comm_object new
  void * operator new(size_t s, int i) {return comm_object::operator new(s,i);}
  void * operator new(size_t s) { return comm_object::operator new(s); }
  void * operator new(size_t, void *ptr) { return ptr; }

  // pack and unpack functions
  void * pack (int *length);
  void unpack (void *in);
};


class ComputeDPMEResultsMsg : public comm_object {
public:
  // data members
  int node;
  int numParticles;
  PmeVector *forces;

  // constructor and destructor
  ComputeDPMEResultsMsg(void);
  ~ComputeDPMEResultsMsg(void);

  // standard new overload for comm_object new
  void * operator new(size_t s, int i) {return comm_object::operator new(s,i);}
  void * operator new(size_t s) { return comm_object::operator new(s); }
  void * operator new(size_t, void *ptr) { return ptr; }

  // pack and unpack functions
  void * pack (int *length);
  void unpack (void *in);
};


#endif


/***************************************************************************
 * RCS INFORMATION:
 *
 *	$RCSfile $
 *	$Author $	$Locker:  $		$State: Exp $
 *	$Revision: 1.1 $	$Date: 1998/04/10 04:15:58 $
 *
 ***************************************************************************
 * REVISION HISTORY:
 *
 * $Log: ComputeDPMEMsgs.h,v $
 * Revision 1.1  1998/04/10 04:15:58  jim
 * Finished incorporating DPME.
 *
 *
 *
 ***************************************************************************/
