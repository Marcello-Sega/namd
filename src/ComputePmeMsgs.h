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
  PmeVector *forces;

  // constructor and destructor
  ComputePmeResultsMsg(void);
  ~ComputePmeResultsMsg(void);

  // pack and unpack functions
  static void* pack(ComputePmeResultsMsg *msg);
  static ComputePmeResultsMsg* unpack(void *ptr);
};


#endif


/***************************************************************************
 * RCS INFORMATION:
 *
 *	$RCSfile $
 *	$Author $	$Locker:  $		$State: Exp $
 *	$Revision: 1.1 $	$Date: 1999/06/08 14:52:07 $
 *
 ***************************************************************************
 * REVISION HISTORY:
 *
 * $Log: ComputePmeMsgs.h,v $
 * Revision 1.1  1999/06/08 14:52:07  jim
 * Incorporated Justin's faster PME code along side DPME.
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
