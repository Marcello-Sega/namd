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

#include "ckdefs.h"
#include "chare.h"
#include "c++interface.h"

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



/***************************************************************************
 * RCS INFORMATION:
 *
 *	$RCSfile $
 *	$Author $	$Locker:  $		$State: Exp $
 *	$Revision: 1.1 $	$Date: 1997/03/19 11:53:51 $
 *
 ***************************************************************************
 * REVISION HISTORY:
 *
 * $Log: BroadcastClient.h,v $
 * Revision 1.1  1997/03/19 11:53:51  ari
 * Add Broadcast mechanism.
 * Fixed RCS Log entries on files that did not have Log entries.
 * Added some register variables to Molecule and ComputeNonbondedExcl.C
 *
 ***************************************************************************/
