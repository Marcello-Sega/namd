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
  ComputePatchPair(ComputeID c, PatchID pid[], int t[]);
  virtual ~ComputePatchPair();

  virtual void initialize();
  virtual void atomUpdate();
  virtual void doWork();

protected :
  int numAtoms[2];
  virtual void doForce(Position* p[2], Force* f[2], AtomProperties* a[2]);

private:
  PatchID patchID[2];
  int trans[2];
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
 *	$Revision: 1.1000 $	$Date: 1997/02/06 15:58:18 $
 *
 ***************************************************************************
 * REVISION HISTORY:
 *
 * $Log: ComputePatchPair.h,v $
 * Revision 1.1000  1997/02/06 15:58:18  ari
 * Resetting CVS to merge branches back into the main trunk.
 * We will stick to main trunk development as suggested by CVS manual.
 * We will set up tags to track fixed points of development/release
 * as suggested by CVS manual - all praise the CVS manual.
 *
 * Revision 1.779  1997/02/06 15:53:10  ari
 * Updating Revision Line, getting rid of branches
 *
 * Revision 1.778.2.2  1997/02/06 02:35:22  jim
 * Implemented periodic boundary conditions - may not work with
 * atom migration yet, but doesn't seem to alter calculation,
 * appears to work correctly when turned on.
 * NamdState chdir's to same directory as config file in argument.
 *
 * Revision 1.778.2.1  1997/02/05 22:18:12  ari
 * Added migration code - Currently the framework is
 * there with compiling code.  This version does
 * crash shortly after migration is complete.
 * Migration appears to complete, but Patches do
 * not appear to be left in a correct state.
 *
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

