//-*-c++-*-
/***************************************************************************/
/*                                                                         */
/*              (C) Copyright 1996 The Board of Trustees of the            */
/*                          University of Illinois                         */
/*                           All Rights Reserved                           */
/*									   */
/***************************************************************************/

/***************************************************************************
 * DESCRIPTION:
 *
 ***************************************************************************/

#ifndef COMPUTEPPAIR_H
#define COMPUTEPPAIR_H

#include "Compute.h"

#include "Templates/Box.h"
#include "Templates/OwnerBox.h"
#include "PositionBox.h"
#include "PositionOwnerBox.h"

class Patch;
class Node;
class PatchMap;

class ComputePatchPair : public Compute {

public:
  ComputePatchPair(ComputeID c, PatchID pid[]);
  virtual ~ComputePatchPair();

  virtual void mapReady();
  virtual void doWork();

protected :
  int numAtoms[2];
  virtual void doForce(Position* p[2], Force* f[2], AtomProperties* a[2]);

private:
  PatchID patchID[2];
  Patch *patch[2];
  PositionBox<Patch> *positionBox[2];
  Box<Patch,Force> *forceBox[2];
  Box<Patch,AtomProperties> *atomBox[2];

};

#endif
/***************************************************************************
 * RCS INFORMATION:
 *
 *	$RCSfile: ComputePatchPair.h,v $
 *	$Author: ari $	$Locker:  $		$State: Exp $
 *	$Revision: 1.778 $	$Date: 1997/01/28 00:30:29 $
 *
 ***************************************************************************
 * REVISION HISTORY:
 *
 * $Log: ComputePatchPair.h,v $
 * Revision 1.778  1997/01/28 00:30:29  ari
 * internal release uplevel to 1.778
 *
 * Revision 1.777.2.2  1997/01/27 22:45:10  ari
 * Basic Atom Migration Code added.
 * Added correct magic first line to .h files for xemacs to go to C++ mode.
 * Compiles and runs without migration turned on.
 *
 * Revision 1.777.2.1  1997/01/24 22:00:31  jim
 * Changes for periodic boundary conditions.
 *
 * Revision 1.777  1997/01/17 19:36:04  ari
 * Internal CVS leveling release.  Start development code work
 * at 1.777.1.1.
 *
 * Revision 1.5  1996/11/23 22:59:57  jim
 * made mapReady() public
 *
 * Revision 1.4  1996/10/30 01:16:32  jim
 * added AtomProperties structure in Patch plus boxes, passing, etc.
 *
 * Revision 1.3  1996/10/30 00:16:16  jim
 * Removed PositionArray usage.
 *
 * Revision 1.2  1996/10/29 23:53:58  jim
 * cleaned up, now only compile blocks are PatchMap, Patch, Compute.
 *
 * Revision 1.1  1996/10/29 22:43:35  ari
 * Initial revision
 *
 * Revision 1.3  1996/10/16 08:22:39  ari
 * *** empty log message ***
 *
 * Revision 1.1  1996/08/19 22:07:49  ari
 * Initial revision
 *
 * Revision 1.4  1996/07/16 01:54:12  ari
 * *** empty log message ***
 *
 * Revision 1.3  96/07/16  01:10:26  01:10:26  ari (Aritomo Shinozaki)
 * Fixed comments, added methods
 * 
 * Revision 1.2  1996/06/25 21:10:48  gursoy
 * *** empty log message ***
 *
 * Revision 1.1  1996/06/24 14:12:26  gursoy
 * Initial revision
 *
 ***************************************************************************/

