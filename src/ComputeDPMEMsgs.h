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


/***************************************************************************
 * RCS INFORMATION:
 *
 *	$RCSfile $
 *	$Author $	$Locker:  $		$State: Exp $
 *	$Revision: 1.3 $	$Date: 1999/09/08 16:05:44 $
 *
 ***************************************************************************
 * REVISION HISTORY:
 *
 * $Log: ComputeDPMEMsgs.h,v $
 * Revision 1.3  1999/09/08 16:05:44  jim
 * Added internal PUB3DFFT package.
 *
 * Revision 1.2  1999/05/11 23:56:18  brunner
 * Changes for new charm version
 *
 * Revision 1.1  1998/04/10 04:15:58  jim
 * Finished incorporating DPME.
 *
 *
 *
 ***************************************************************************/
