//-*-c++-*-
/***************************************************************************/
/*         (C) Copyright 1996,1997 The Board of Trustees of the            */
/*                          University of Illinois                         */
/*                           All Rights Reserved                           */
/***************************************************************************/
/***************************************************************************
 * DESCRIPTION:  Messages needed for sending SMD data.
 *
 ***************************************************************************/

#ifndef SMDDATAMSGS_H
#define SMDDATAMSGS_H

#include "charm++.h"

#include "NamdTypes.h"
#include "Node.decl.h"

class SMDDataMsg : public CMessage_SMDDataMsg {
public:
  // data members
  int curTime;            // current timestep
  int timeStamp;          // time the data was last updated
  Vector direction;       // direction of restraint point movement
  Vector refPos;          // restraint point position
  Vector atomPosVmin;     // restrained atom position for Vmin averaging
  Vector atomPosVmax;     // restrained atom position for Vmax averaging

  // constructor and destructor
  SMDDataMsg(void) {};
  ~SMDDataMsg(void) {}; 

  // pack and unpack functions
  static void* pack(SMDDataMsg *msg);
  static SMDDataMsg* unpack(void *ptr);
};


#endif // SMDDATAMSGS_H



