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

#ifndef HOMEPATCHTYPES_H
#define HOMEPATCHTYPES_H

#include "PatchTypes.h"
#include "Templates/Box.h"

class ProxyListElem {
public:
  ProxyListElem() {};
  ProxyListElem(NodeID n) : node(n), forceBox(0) {};
  ProxyListElem(NodeID n, Box<Patch,Results> *f ) : node(n), forceBox(f) {};
  ~ProxyListElem() {};

  int operator==(const ProxyListElem &p) {
    return (node == p.node);
  }

  NodeID node;
  Box<Patch,Results> *forceBox;
};

typedef ResizeArray<ProxyListElem> ProxyList;
typedef ResizeArrayIter<ProxyListElem> ProxyListIter;


#endif
/***************************************************************************
 * RCS INFORMATION:
 *
 *	$RCSfile $
 *	$Author $	$Locker:  $		$State: Exp $
 *	$Revision: 1.1003 $	$Date: 1997/07/09 21:26:41 $
 *
 ***************************************************************************
 * REVISION HISTORY:
 *
 * $Log: HomePatchTypes.h,v $
 * Revision 1.1003  1997/07/09 21:26:41  milind
 * Ported NAMD2 to SP3. The SP specific code is within #ifdef SP2
 * and #endif's.
 *
 * Revision 1.1002  1997/03/19 11:54:18  ari
 * Add Broadcast mechanism.
 * Fixed RCS Log entries on files that did not have Log entries.
 * Added some register variables to Molecule and ComputeNonbondedExcl.C
 *
 ***************************************************************************/
