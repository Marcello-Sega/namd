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
#include "BroadcastObject.h"

template <class T>
T
SimpleBroadcastObject<T>::get(int tag) {
  void *buf;

  while (!(buf = (BroadcastMgr::Object())->getbuf(*this, tag))) {
    suspendFor(tag);
  }
  T tmp = *(T *)buf;
  delete buf;
  return tmp;
}
    
template <class T>
void 
SimpleBroadcastObject<T>::publish(int tag,const T &t ) {
  void *buf = new T;
  *(T *)buf = t;
  BroadcastMgr::Object()->send(*this, tag, buf, sizeof(T));
}

/***************************************************************************
 * RCS INFORMATION:
 *
 *	$RCSfile $
 *	$Author $	$Locker:  $		$State: Exp $
 *	$Revision: 1.1 $	$Date: 1997/03/19 11:53:54 $
 *
 ***************************************************************************
 * REVISION HISTORY:
 *
 * $Log: BroadcastObject.C,v $
 * Revision 1.1  1997/03/19 11:53:54  ari
 * Add Broadcast mechanism.
 * Fixed RCS Log entries on files that did not have Log entries.
 * Added some register variables to Molecule and ComputeNonbondedExcl.C
 *
 *
 ***************************************************************************/
