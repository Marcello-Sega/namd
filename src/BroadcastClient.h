//-*-c++-*-
/***************************************************************************/
/*              (C) Copyright 1996,1997 The Board of Trustees of the       */
/*                          University of Illinois                         */
/*                           All Rights Reserved                           */
/***************************************************************************/
/***************************************************************************
 * DESCRIPTION: Coordinates broadcast of a data type from a Controller/Seq
 *		to all other Controller/Sequencer type objects (they must
 *		run in a thread!)
 ***************************************************************************/

#include "charm++.h"

#ifndef _BCASTCLI_H
#define _BCASTCLI_H

class BroadcastClient {
public:
  int id;
  BroadcastClient(int id);
  ~BroadcastClient();
  void awaken(int id, int tag);

protected: 
  void suspendFor(int tag);

  int suspended;
  int waitForTag;
  CthThread thread;
};

#endif

