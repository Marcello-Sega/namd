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
#include "BroadcastMgr.h"
#include "BroadcastClient.h"

#ifndef _BCASTOBJ_H
#define _BCASTOBJ_H

template <class T>
class SimpleBroadcastObject : public BroadcastClient {
public:
  SimpleBroadcastObject(int id=0) : BroadcastClient(id) { }
  ~SimpleBroadcastObject() { }

  T get(int);
  void publish(int, const T &);
};

#endif


/***************************************************************************
 * RCS INFORMATION:
 *
 *	$RCSfile $
 *	$Author $	$Locker:  $		$State: Exp $
 *	$Revision: 1.1 $	$Date: 1997/03/19 11:53:55 $
 *
 ***************************************************************************
 * REVISION HISTORY:
 *
 * $Log: BroadcastObject.h,v $
 * Revision 1.1  1997/03/19 11:53:55  ari
 * Add Broadcast mechanism.
 * Fixed RCS Log entries on files that did not have Log entries.
 * Added some register variables to Molecule and ComputeNonbondedExcl.C
 *
 *
 ***************************************************************************/
