/***************************************************************************/
/*     (C) Copyright 1996,1997 The Board of Trustees of the                */
/*                          University of Illinois                         */
/*                           All Rights Reserved                           */
/***************************************************************************/
/***************************************************************************
 * DESCRIPTION: Top of Compute hierarchy.  
 *	enqueueWork() - delivers Compute object itself to queue up for
 *			doWork()
 *	doWork() - called by work queue
 ***************************************************************************/

#include "main.h"
#include "ckdefs.h"
#include "chare.h"
#include "c++interface.h"

#include "WorkDistrib.top.h"
#include "WorkDistrib.h"

#include "NamdTypes.h"
#include "Templates/Box.h"
#include "Templates/OwnerBox.h"

#include "Node.h"
#include "Compute.h"

#define MIN_DEBUG_LEVEL 4
#define DEBUGM
#include "Debug.h"

Node *Compute::node=0;
PatchMap *Compute::patchMap=0;

int Compute::totalComputes = 0;

Compute::Compute(ComputeID c) : cid(c), basePriority(DEFPRIO){ 
  totalComputes++;
  doAtomUpdate = false;
}

void Compute::enqueueWork() {
  if (!this) { iout << iPE << iERRORF << "This Compute is NULL!!!\n" << endi; }
  if ( ! noWork() )
  {
    WorkDistrib::messageEnqueueWork(this);  // should be in ComputeMgr?
  }
}

//---------------------------------------------------------------------
// Signal from patch or proxy that data is ready.
// When all Patches and Proxies needed by this Compute object
// have checked-in, we are ready to enqueueWork()
//---------------------------------------------------------------------
void Compute::patchReady(PatchID patchID, int doneMigration) { 
  if (doneMigration) { // If any patch has done migration - we must remap
    doAtomUpdate = true; 
  }

  if (numPatches <= 0) {
    iout << iERRORF 
      << "Compute::patchReady("<<patchID<<")-call not valid!\n"
      << endi;
  } else {
    if (! --patchReadyCounter) {
      patchReadyCounter = numPatches;
      if (doAtomUpdate) {
	atomUpdate();
	doAtomUpdate = false;
      }
      enqueueWork();
    }
  }
}


int Compute::noWork() {
  return 0;
}

void Compute::doWork() {
  iout << iERRORF 
    << "Default Compute::doWork() called.\n"
    << endi;
}

int Compute::sequence(void)
{
  return -1;
}

/***************************************************************************
 * RCS INFORMATION:
 *
 *	$RCSfile: Compute.C,v $
 *	$Author: milind $	$Locker:  $		$State: Exp $
 *	$Revision: 1.1014 $	$Date: 1997/09/28 10:19:03 $
 *
 ***************************************************************************
 * REVISION HISTORY:
 *
 * $Log: Compute.C,v $
 * Revision 1.1014  1997/09/28 10:19:03  milind
 * Fixed priorities, ReductionMgr etc.
 *
 * Revision 1.1013  1997/08/26 16:26:11  jim
 * Revamped prioritites for petter performance and easier changes.
 *
 * Revision 1.1012  1997/08/20 23:27:36  jim
 * Created multiple enqueueWork entry points to aid analysis.
 *
 * Revision 1.1011  1997/04/10 09:13:47  ari
 * Final debugging for compute migration / proxy creation for load balancing.
 * Lots of debug code added, mostly turned off now.
 * Fixed bug in PositionBox when Patch had no dependencies.
 * Eliminated use of cout and misuse of iout in numerous places.
 *                                            Ari & Jim
 *
 * Revision 1.1010  1997/04/08 07:08:08  ari
 * Modification for dynamic loadbalancing - moving computes
 * Still bug in new computes or usage of proxies/homepatches.
 * Works if ldbStrategy is none as before.
 *
 * Revision 1.1009  1997/04/06 22:44:57  ari
 * Add priorities to messages.  Mods to help proxies without computes.
 * Added quick enhancement to end of list insertion of ResizeArray(s)
 *
 * Revision 1.1008  1997/04/04 23:34:14  milind
 * Got NAMD2 to run on Origin2000.
 * Included definitions of class static variables in C files.
 * Fixed alignment bugs by using memcpy instead of assignment in
 * pack and unpack.
 *
 * Revision 1.1007  1997/04/03 23:22:16  jim
 * Added basic priority() method to Compute.  Only distinguishes between
 * local and nonlocal computations for now.
 *
 * Revision 1.1006  1997/03/20 23:53:26  ari
 * Some changes for comments. Copyright date additions.
 * Hooks for base level update of Compute objects from ComputeMap
 * by ComputeMgr.  Useful for new compute migration functionality.
 *
 * Revision 1.1005  1997/03/19 05:49:46  jim
 * Added ComputeSphericalBC, cleaned up make dependencies.
 *
 * Revision 1.1004  1997/03/12 23:59:37  jim
 * Added Compute::noWork() protocol to not enqueue do-nothing compute objects.
 *
 * Revision 1.1003  1997/03/06 22:05:57  ari
 * Removed Compute.ci
 * Comments added - more code cleaning
 *
 * Revision 1.1002  1997/02/13 16:17:11  ari
 * Intermediate debuging commit - working to fix deep bug in migration?
 *
 * Revision 1.1001  1997/02/11 18:51:39  ari
 * Modified with #ifdef DPMTA to safely eliminate DPMTA codes
 * fixed non-buffering of migration msgs
 * Migration works on multiple processors
 *
 * Revision 1.1000  1997/02/06 15:57:41  ari
 * Resetting CVS to merge branches back into the main trunk.
 * We will stick to main trunk development as suggested by CVS manual.
 * We will set up tags to track fixed points of development/release
 * as suggested by CVS manual - all praise the CVS manual.
 *
 * Revision 1.779  1997/02/06 15:52:48  ari
 * Updating Revision Line, getting rid of branches
 *
 * Revision 1.778.2.1  1997/02/05 22:17:54  ari
 * Added migration code - Currently the framework is
 * there with compiling code.  This version does
 * crash shortly after migration is complete.
 * Migration appears to complete, but Patches do
 * not appear to be left in a correct state.
 *
 * Revision 1.778  1997/01/28 00:29:59  ari
 * internal release uplevel to 1.778
 *
 * Revision 1.777.2.1  1997/01/27 22:44:55  ari
 * Basic Atom Migration Code added.
 * Added correct magic first line to .h files for xemacs to go to C++ mode.
 * Compiles and runs without migration turned on.
 *
 * Revision 1.777  1997/01/17 19:35:32  ari
 * Internal CVS leveling release.  Start development code work
 * at 1.777.1.1.
 *
 * Revision 1.4  1996/11/30 20:30:36  jim
 * turned off some debugging, switched to DebugM()
 *
 * Revision 1.3  1996/11/22 00:18:51  ari
 * *** empty log message ***
 *
 * Revision 1.2  1996/10/22 19:12:16  ari
 * *** empty log message ***
 *
 * Revision 1.1  1996/10/16 08:22:39  ari
 * Initial revision
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

