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


static char ident[] = "@(#)$Header: /home/cvs/namd/cvsroot/namd2/src/ProxyPatch.C,v 1.1018 1998/03/03 23:05:25 brunner Exp $";

#include "charm++.h"

#include "main.top.h"
#include "main.h"
#include "ProxyPatch.h"
#include "ProxyMgr.top.h"
#include "ProxyMgr.h"
#include "AtomMap.h"

#define MIN_DEBUG_LEVEL 4
//#define  DEBUGM
#include "Debug.h"

ProxyPatch::ProxyPatch(PatchID pd) : 
  Patch(pd), msgBuffer(NULL), msgAllBuffer(NULL)
{
  DebugM(4, "ProxyPatch(" << pd << ") at " << this << "\n");
  ProxyMgr::Object()->registerProxy(patchID);
}

void ProxyPatch::boxClosed(int box)
{
  if ( box == 1 ) {
    sendResults();
  }
  if ( ! --boxesOpen ) {
    DebugM(2,patchID << ": " << "Checking message buffer.\n");
    if ( msgBuffer ) {
      DebugM(3,"Patch " << patchID << " processing buffered proxy data.\n");
      receiveData(msgBuffer);
    } else if (msgAllBuffer ) {
      DebugM(3,"Patch " << patchID << " processing buffered proxy ALL data.\n");
      receiveAll(msgAllBuffer);
    }
  }
  else {
    DebugM(3,"ProxyPatch " << patchID << ": " << boxesOpen << " boxes left to close.\n");
  }
}

void ProxyPatch::receiveAtoms(ProxyAtomsMsg *msg)
{
  DebugM(3, "receiveAtoms(" << patchID << ")\n");
  loadAtoms(msg->atomIDList);
  AtomMap::Object()->registerIDs(patchID,msg->atomIDList);
  delete msg;
}

void ProxyPatch::receiveData(ProxyDataMsg *msg)
{
  DebugM(3, "receiveData(" << patchID << ")\n");
  if ( boxesOpen )
  {
    // store message in queue (only need one element, though)
    msgBuffer = msg;
    return;
  }
  msgBuffer = NULL;
  flags = msg->flags;
  p = msg->positionList;
  delete msg;
  positionsReady(0);
}

void ProxyPatch::receiveAll(ProxyAllMsg *msg)
{
  DebugM(3, "receiveData(" << patchID << ")\n");
  if ( boxesOpen )
  {
    // store message in queue (only need one element, though)
    msgAllBuffer = msg;
    return;
  }
  msgAllBuffer = NULL;

  AtomMap::Object()->unregisterIDs(patchID,atomIDList);
  loadAtoms(msg->atomIDList);
  AtomMap::Object()->registerIDs(patchID,msg->atomIDList);
  flags = msg->flags;
  p = msg->positionList;

  delete msg;

  positionsReady(1);
}

void ProxyPatch::sendResults(void)
{
  DebugM(3, "sendResults(" << patchID << ")\n");
  ProxyResultMsg *msg 
    = new (MsgIndex(ProxyResultMsg)) ProxyResultMsg;
  msg->node = CMyPe();
  msg->patch = patchID;
  register int i = 0;
  register ForceList::iterator f_i, f_e, f2_i;
  for ( i = Results::normal + 1 ; i <= flags.maxForceMerged; ++i ) {
    f_i = f[Results::normal].begin(); f_e = f[Results::normal].end();
    f2_i = f[i].begin();
    for ( ; f_i != f_e; ++f_i, ++f2_i ) *f_i += *f2_i;
    f[i].resize(0);
  }
  for ( i = flags.maxForceUsed + 1; i < Results::maxNumForces; ++i )
    f[i].resize(0);
  for ( i = 0; i < Results::maxNumForces; ++i ) 
    msg->forceList[i] = f[i];
  ProxyMgr::Object()->sendResults(msg);
}

/***************************************************************************
 * RCS INFORMATION:
 *
 *	$RCSfile: ProxyPatch.C,v $
 *	$Author: brunner $	$Locker:  $		$State: Exp $
 *	$Revision: 1.1018 $	$Date: 1998/03/03 23:05:25 $
 *
 ***************************************************************************
 * REVISION HISTORY:
 *
 * $Log: ProxyPatch.C,v $
 * Revision 1.1018  1998/03/03 23:05:25  brunner
 * Changed include files for new simplified Charm++ include file structure.
 *
 * Revision 1.1017  1997/12/22 21:29:27  jim
 * Proxies no longer send empty arrays back to HomePatch.  Requires some new
 * flags to be set correctly in Sequencer in order to work.  These are:
 *   maxForceMerged - this and faster are added into Results::normal array
 *   maxForceUsed - all forces slower than this are discarded (assumed zero)
 * Generally maxForceMerged doesn't change but maxForceUsed depends on timestep.
 *
 * Revision 1.1016  1997/09/28 10:19:07  milind
 * Fixed priorities, ReductionMgr etc.
 *
 * Revision 1.1015  1997/08/26 16:26:16  jim
 * Revamped prioritites for petter performance and easier changes.
 *
 * Revision 1.1014  1997/08/22 20:12:04  milind
 * Turned on Priorities.
 *
 * Revision 1.1013  1997/07/08 15:48:11  milind
 * Made namd2 to work with Origin2000: Again...
 *
 * Revision 1.1012  1997/04/10 09:14:10  ari
 * Final debugging for compute migration / proxy creation for load balancing.
 * Lots of debug code added, mostly turned off now.
 * Fixed bug in PositionBox when Patch had no dependencies.
 * Eliminated use of cout and misuse of iout in numerous places.
 *                                            Ari & Jim
 *
 * Revision 1.1011  1997/04/08 07:08:58  ari
 * Modification for dynamic loadbalancing - moving computes
 * Still bug in new computes or usage of proxies/homepatches.
 * Works if ldbStrategy is none as before.
 *
 * Revision 1.1010  1997/04/06 22:45:12  ari
 * Add priorities to messages.  Mods to help proxies without computes.
 * Added quick enhancement to end of list insertion of ResizeArray(s)
 *
 * Revision 1.1009  1997/03/12 22:06:48  jim
 * First step towards multiple force returns and multiple time stepping.
 *
 * Revision 1.1008  1997/02/28 04:47:12  jim
 * Full electrostatics now works with fulldirect on one node.
 *
 * Revision 1.1007  1997/02/26 16:53:17  ari
 * Cleaning and debuging for memory leaks.
 * Adding comments.
 * Removed some dead code due to use of Quiescense detection.
 *
 * Revision 1.1006  1997/02/13 23:17:20  ari
 * Fixed a final bug in AtomMigration - numatoms in ComputePatchPair.C not
 * set correctly in atomUpdate()
 *
 * Revision 1.1005  1997/02/13 16:17:19  ari
 * Intermediate debuging commit - working to fix deep bug in migration?
 *
 * Revision 1.1004  1997/02/10 08:26:03  jim
 * Turned off debugging.
 *
 * Revision 1.1003  1997/02/07 17:39:41  ari
 * More debugging for atomMigration.
 * Using -w on CC got us some minor fixes
 * using purify got us a major memory problem due to bad sizing of dummy force
 *
 * Revision 1.1002  1997/02/07 16:49:29  jim
 * Fixing bugs that affect parallel atom migration.
 *
 * Revision 1.1001  1997/02/07 05:42:32  ari
 * Some bug fixing - atom migration on one node works
 * Atom migration on multiple nodes gets SIGSEGV
 *
 * Revision 1.1000  1997/02/06 15:59:13  ari
 * Resetting CVS to merge branches back into the main trunk.
 * We will stick to main trunk development as suggested by CVS manual.
 * We will set up tags to track fixed points of development/release
 * as suggested by CVS manual - all praise the CVS manual.
 *
 * Revision 1.779  1997/02/06 15:53:27  ari
 * Updating Revision Line, getting rid of branches
 *
 * Revision 1.778.2.1  1997/02/05 22:18:22  ari
 * Added migration code - Currently the framework is
 * there with compiling code.  This version does
 * crash shortly after migration is complete.
 * Migration appears to complete, but Patches do
 * not appear to be left in a correct state.
 *
 * Revision 1.778  1997/01/28 00:31:19  ari
 * internal release uplevel to 1.778
 *
 * Revision 1.777.2.2  1997/01/27 22:45:39  ari
 * Basic Atom Migration Code added.
 * Added correct magic first line to .h files for xemacs to go to C++ mode.
 * Compiles and runs without migration turned on.
 *
 * Revision 1.777.2.1  1997/01/22 20:25:52  nealk
 * Disabled debugging.
 *
 * Revision 1.777  1997/01/17 19:36:54  ari
 * Internal CVS leveling release.  Start development code work
 * at 1.777.1.1.
 *
 * Revision 1.11  1996/12/19 00:33:10  jim
 * added message buffering debugging
 *
 * Revision 1.10  1996/12/17 23:58:02  jim
 * proxy result reporting is working
 *
 * Revision 1.9  1996/12/17 22:13:22  jim
 * implemented ProxyDataMsg use
 *
 * Revision 1.8  1996/12/16 22:52:43  jim
 * added placement new and explicit destructor calls to ProxyAtomsMsg
 *
 * Revision 1.7  1996/12/14 00:02:42  jim
 * debugging ProxyAtomsMsg path to make compute creation work
 *
 * Revision 1.6  1996/12/05 23:45:09  ari
 * *** empty log message ***
 *
 * Revision 1.5  1996/12/05 22:17:48  jim
 * removed .bot.h file
 *
 * Revision 1.4  1996/12/05 22:09:44  jim
 * fixed compile errors
 *
 * Revision 1.3  1996/12/05 22:02:17  jim
 * added positionsReady call to receiveData
 *
 * Revision 1.2  1996/12/05 21:11:06  jim
 * filled out message functions
 *
 * Revision 1.1  1996/12/05 01:44:16  ari
 * Initial revision
 *
 * Revision 1.13  1996/12/01 02:31:37  jim
 * improved debugging, fixed boxesOpen possible bug
 *
 * Revision 1.12  1996/11/30 00:41:24  jim
 * added boxesOpen counting to support HomePatch::boxClosed()
 *
 * Revision 1.11  1996/11/22 01:02:18  ari
 * *** empty log message ***
 *
 * Revision 1.10  1996/11/22 00:18:51  ari
 * *** empty log message ***
 *
 * Revision 1.9  1996/11/04 17:13:31  ari
 * *** empty log message ***
 *
 * Revision 1.8  1996/11/01 21:20:45  ari
 * *** empty log message ***
 *
 * Revision 1.7  1996/10/30 01:50:19  jim
 * atom properties list now filled on creation
 *
 * Revision 1.6  1996/10/30 01:16:32  jim
 * added AtomProperties structure in Patch plus boxes, passing, etc.
 *
 * Revision 1.5  1996/10/29 23:35:27  ari
 * *** empty log message ***
 *
 * Revision 1.4  1996/10/16 08:22:39  ari
 * *** empty log message ***
 *
 * Revision 1.3  1996/10/04 21:07:46  jim
 * Moved in functionality from HomePatch
 *
 * Revision 1.2  1996/09/10 03:07:04  ari
 * *** empty log message ***
 *
 * Revision 1.1  1996/08/19 22:07:49  ari
 * Initial revision
 *
 * Revision 1.6  1996/07/15 22:56:48  gursoy
 * *** empty log message ***
 *
 * Revision 1.5  1996/07/15 21:17:55  gursoy
 * *** empty log message ***
 *
 * Revision 1.4  1996/07/09 21:54:34  gursoy
 * *** empty log message ***
 *
 * Revision 1.3  1996/06/24 14:14:30  gursoy
 * *** empty log message ***
 *
 * Revision 1.2  1996/06/12 16:34:23  brunner
 * *** empty log message ***
 *
 * Revision 1.1  1996/05/30 21:31:36  gursoy
 * Initial revision
 *
 ***************************************************************************/
