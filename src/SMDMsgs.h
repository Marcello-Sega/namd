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

class SMDDataMsg : public comm_object {
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

  // standard new overload for comm_object new
  void * operator new(size_t s, int i) {return comm_object::operator new(s,i);}
  void * operator new(size_t s) { return comm_object::operator new(s); }
  void * operator new(size_t, void *ptr) { return ptr; }

  // pack and unpack functions
  void * pack (int *length);
  void unpack (void *in);
};


#endif // SMDDATAMSGS_H



