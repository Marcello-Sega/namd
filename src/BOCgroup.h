//-*-c++-*-
/***************************************************************************/
/*                                                                         */
/*              (C) Copyright 1996 The Board of Trustees of the            */
/*                          University of Illinois                         */
/*                           All Rights Reserved                           */
/*                                                                         */
/***************************************************************************/

/***************************************************************************
 * DESCRIPTION:
 *
 ***************************************************************************/

#ifndef BOCGROUP_H
#define BOCGROUP_H

#include "charm++.h"

class BOCgroup {
public:
  int workDistrib;
  int patchMgr;
  int proxyMgr;
  int computeMgr;
  int reductionMgr;
  int collectionMgr;
  int broadcastMgr;
  int ldbCoordinator;
  int node;
};

class BOCclass : public Group {
};

#endif /* BOCGROUP_H */


/***************************************************************************
 * RCS INFORMATION:
 *
 *	$RCSfile $
 *	$Author $	$Locker:  $		$State: Exp $
 *	$Revision: 1.1005 $	$Date: 1999/05/11 23:56:12 $
 *
 ***************************************************************************
 * REVISION HISTORY:
 *
 * $Log: BOCgroup.h,v $
 * Revision 1.1005  1999/05/11 23:56:12  brunner
 * Changes for new charm version
 *
 * Revision 1.1004  1998/03/03 23:05:00  brunner
 * Changed include files for new simplified Charm++ include file structure.
 *
 * Revision 1.1003  1997/11/07 20:17:31  milind
 * Made NAMD to run on shared memory machines.
 *
 * Revision 1.1002  1997/03/27 20:25:39  brunner
 * Changes for LdbCoordinator, the load balance control BOC
 *
 * Revision 1.1001  1997/03/19 11:53:50  ari
 * Add Broadcast mechanism.
 * Fixed RCS Log entries on files that did not have Log entries.
 * Added some register variables to Molecule and ComputeNonbondedExcl.C
 *
 ***************************************************************************/
