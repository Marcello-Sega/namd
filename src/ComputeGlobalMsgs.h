//-*-c++-*-
/***************************************************************************/
/*         (C) Copyright 1996,1997 The Board of Trustees of the            */
/*                          University of Illinois                         */
/*                           All Rights Reserved                           */
/***************************************************************************/
/***************************************************************************
 * DESCRIPTION:  Messages needed for ComputeGlobal operation.
 *
 ***************************************************************************/

#ifndef COMPUTEGLOBALMSGS_H
#define COMPUTEGLOBALMSGS_H

#include "ckdefs.h"
#include "chare.h"
#include "c++interface.h"

#include "NamdTypes.h"

class ComputeGlobalConfigMsg : public comm_object {
public:
  // data members
  AtomIDList aid;

  // constructor and destructor
  ComputeGlobalConfigMsg(void);
  ~ComputeGlobalConfigMsg(void);

  // standard new overload for comm_object new
  void * operator new(size_t s, int i) {return comm_object::operator new(s,i);}
  void * operator new(size_t s) { return comm_object::operator new(s); }
  void * operator new(size_t, void *ptr) { return ptr; }

  // pack and unpack functions
  void * pack (int *length);
  void unpack (void *in);
};


class ComputeGlobalDataMsg : public comm_object {
public:
  // data members
  AtomIDList aid;
  PositionList p;

  // constructor and destructor
  ComputeGlobalDataMsg(void);
  ~ComputeGlobalDataMsg(void);

  // standard new overload for comm_object new
  void * operator new(size_t s, int i) {return comm_object::operator new(s,i);}
  void * operator new(size_t s) { return comm_object::operator new(s); }
  void * operator new(size_t, void *ptr) { return ptr; }

  // pack and unpack functions
  void * pack (int *length);
  void unpack (void *in);
};


class ComputeGlobalResultsMsg : public comm_object {
public:
  // data members
  AtomIDList aid;
  ForceList f;
  int reconfig;
  AtomIDList newaid;

  // constructor and destructor
  ComputeGlobalResultsMsg(void);
  ~ComputeGlobalResultsMsg(void);

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
 *	$Revision: 1.1 $	$Date: 1997/12/19 23:48:47 $
 *
 ***************************************************************************
 * REVISION HISTORY:
 *
 * $Log: ComputeGlobalMsgs.h,v $
 * Revision 1.1  1997/12/19 23:48:47  jim
 * Added Tcl interface for calculating forces.
 *
 *
 ***************************************************************************/
