//-*-c++-*-
/***************************************************************************/
/*                                                                         */
/*              (C) Copyright 1996 The Board of Trustees of the            */
/*                          University of Illinois                         */
/*                           All Rights Reserved                           */
/*									   */
/***************************************************************************/

/***************************************************************************
 * DESCRIPTION: ComputeNonbondedPair.h
 *
 ***************************************************************************/

#ifndef COMPUTENONBONDEDPAIR_H
#define COMPUTENONBONDEDPAIR_H

#include "ComputePatchPair.h"
#include "ComputeNonbondedUtil.h"

class ComputeNonbondedPair : public ComputePatchPair, private ComputeNonbondedUtil {

public:
  ComputeNonbondedPair(ComputeID c, PatchID pid[], int trans[]);
  ~ComputeNonbondedPair();

protected :
  virtual void initialize();
  virtual int noWork();
  virtual void doForce(Position* p[2], Results* r[2], AtomProperties* a[2]);

  PositionBox<Patch> *avgPositionBox[2];

  SubmitReduction *reduction;

};

#endif
/***************************************************************************
 * RCS INFORMATION:
 *
 *	$RCSfile: ComputeNonbondedPair.h,v $
 *	$Author: jim $	$Locker:  $		$State: Exp $
 *	$Revision: 1.1009 $	$Date: 1999/08/20 19:11:10 $
 *
 ***************************************************************************
 * REVISION HISTORY:
 *
 * $Log: ComputeNonbondedPair.h,v $
 * Revision 1.1009  1999/08/20 19:11:10  jim
 * Added MOLLY - mollified impluse method.
 *
 * Revision 1.1008  1999/06/17 15:46:10  jim
 * Completely rewrote reduction system to eliminate need for sequence numbers.
 *
 * Revision 1.1007  1997/08/26 16:26:13  jim
 * Revamped prioritites for petter performance and easier changes.
 *
 * Revision 1.1006  1997/04/08 07:08:23  ari
 * Modification for dynamic loadbalancing - moving computes
 * Still bug in new computes or usage of proxies/homepatches.
 * Works if ldbStrategy is none as before.
 *
 * Revision 1.1005  1997/04/03 23:22:21  jim
 * Added basic priority() method to Compute.  Only distinguishes between
 * local and nonlocal computations for now.
 *
 * Revision 1.1004  1997/03/18 21:35:29  jim
 * Eliminated fake_seq.  Reductions now use Patch::flags.seq.
 *
 * Revision 1.1003  1997/03/13 06:37:01  jim
 * Multiple time-stepping implemented, still needs proper splitting functions.
 *
 * Revision 1.1002  1997/03/12 23:59:42  jim
 * Added Compute::noWork() protocol to not enqueue do-nothing compute objects.
 *
 * Revision 1.1001  1997/02/11 18:51:44  ari
 * Modified with #ifdef DPMTA to safely eliminate DPMTA codes
 * fixed non-buffering of migration msgs
 * Migration works on multiple processors
 *
 * Revision 1.1000  1997/02/06 15:58:11  ari
 * Resetting CVS to merge branches back into the main trunk.
 * We will stick to main trunk development as suggested by CVS manual.
 * We will set up tags to track fixed points of development/release
 * as suggested by CVS manual - all praise the CVS manual.
 *
 * Revision 1.779  1997/02/06 15:53:06  ari
 * Updating Revision Line, getting rid of branches
 *
 * Revision 1.778.2.2  1997/02/06 02:35:20  jim
 * Implemented periodic boundary conditions - may not work with
 * atom migration yet, but doesn't seem to alter calculation,
 * appears to work correctly when turned on.
 * NamdState chdir's to same directory as config file in argument.
 *
 * Revision 1.778.2.1  1997/02/05 22:18:09  ari
 * Added migration code - Currently the framework is
 * there with compiling code.  This version does
 * crash shortly after migration is complete.
 * Migration appears to complete, but Patches do
 * not appear to be left in a correct state.
 *
 * Revision 1.778  1997/01/28 00:30:23  ari
 * internal release uplevel to 1.778
 *
 * Revision 1.777.2.1  1997/01/27 22:45:07  ari
 * Basic Atom Migration Code added.
 * Added correct magic first line to .h files for xemacs to go to C++ mode.
 * Compiles and runs without migration turned on.
 *
 * Revision 1.777  1997/01/17 19:35:57  ari
 * Internal CVS leveling release.  Start development code work
 * at 1.777.1.1.
 *
 * Revision 1.6  1997/01/16 20:00:11  jim
 * Added reduction calls to ComputeNonbondedSelf and ...Pair.
 * Also moved some code from ...Excl to ...Util.
 *
 * Revision 1.5  1996/11/05 22:40:05  jim
 * commented out undefined virtual destructor
 *
 * Revision 1.4  1996/10/31 21:43:29  jim
 * First incarnation as ...Pair
 *
 * Revision 1.3  1996/10/30 01:16:32  jim
 * added AtomProperties structure in Patch plus boxes, passing, etc.
 *
 * Revision 1.2  1996/10/30 00:16:16  jim
 * Removed PositionArray usage.
 *
 * Revision 1.1  1996/10/29 23:55:54  jim
 * Initial revision
 *
 *
 ***************************************************************************/

