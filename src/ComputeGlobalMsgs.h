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

#include "charm++.h"

#include "NamdTypes.h"
#include "ComputeMgr.decl.h"

class ComputeGlobalConfigMsg : public CMessage_ComputeGlobalConfigMsg {
public:
  // data members
  AtomIDList aid;
  AtomIDList gdef;  // group definitions

  // constructor and destructor
  ComputeGlobalConfigMsg(void);
  ~ComputeGlobalConfigMsg(void);

  // pack and unpack functions
  static void* pack(ComputeGlobalConfigMsg *msg);
  static ComputeGlobalConfigMsg* unpack(void *ptr);
};


class ComputeGlobalDataMsg : public CMessage_ComputeGlobalDataMsg {
public:
  // data members
  AtomIDList aid;
  PositionList p;
  PositionList gcom;  // group center of mass

  // constructor and destructor
  ComputeGlobalDataMsg(void);
  ~ComputeGlobalDataMsg(void);

  // pack and unpack functions
  static void* pack(ComputeGlobalDataMsg *msg);
  static ComputeGlobalDataMsg* unpack(void *ptr);
};


class ComputeGlobalResultsMsg : public CMessage_ComputeGlobalResultsMsg {
public:
  // data members
  AtomIDList aid;
  ForceList f;
  ForceList gforce;  // force on group
  int reconfig;
  AtomIDList newaid;
  AtomIDList newgdef;

  // constructor and destructor
  ComputeGlobalResultsMsg(void);
  ~ComputeGlobalResultsMsg(void);

  // pack and unpack functions
  static void* pack(ComputeGlobalResultsMsg *msg);
  static ComputeGlobalResultsMsg* unpack(void *ptr);
};


#endif


/***************************************************************************
 * RCS INFORMATION:
 *
 *	$RCSfile $
 *	$Author $	$Locker:  $		$State: Exp $
 *	$Revision: 1.4 $	$Date: 1999/05/11 23:56:23 $
 *
 ***************************************************************************
 * REVISION HISTORY:
 *
 * $Log: ComputeGlobalMsgs.h,v $
 * Revision 1.4  1999/05/11 23:56:23  brunner
 * Changes for new charm version
 *
 * Revision 1.3  1998/03/03 23:13:47  brunner
 * Changing include files for new charm++ includes
 *
 * Revision 1.2  1998/02/16 00:24:38  jim
 * Added atom group centers of mass to Tcl interface.
 *
 * Revision 1.1  1997/12/19 23:48:47  jim
 * Added Tcl interface for calculating forces.
 *
 *
 ***************************************************************************/
