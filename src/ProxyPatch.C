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


static char ident[] = "@(#)$Header: /home/cvs/namd/cvsroot/namd2/src/ProxyPatch.C,v 1.1001 1997/02/07 05:42:32 ari Exp $";

#include "ckdefs.h"
#include "chare.h"
#include "c++interface.h"

#include "main.top.h"
#include "main.h"
#include "ProxyPatch.h"
#include "ProxyMgr.top.h"
#include "ProxyMgr.h"
#include "AtomMap.h"

#define MIN_DEBUG_LEVEL 4
#define  DEBUGM
#include "Debug.h"

ProxyPatch::ProxyPatch(PatchID pd) : Patch(pd), msgBuffer(NULL)
{
  msgBuffer = NULL;
  msgAllBuffer = NULL;
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
      DebugM(4,"Patch " << patchID << " processing buffered proxy data.\n");
      receiveData(msgBuffer);
    } else if (msgAllBuffer ) {
      DebugM(4,"Patch " << patchID << " processing buffered proxy ALL data.\n");
      receiveAll(msgAllBuffer);
    }
  }
  else {
    DebugM(2,patchID << ": " << boxesOpen << " boxes left to close.\n");
  }
}

void ProxyPatch::receiveAtoms(ProxyAtomsMsg *msg)
{
  loadAtoms(msg->atomIDList);
  loadAtomProperties();
  AtomMap::Object()->registerIDs(patchID,msg->atomIDList);
  delete msg;
}

void ProxyPatch::receiveData(ProxyDataMsg *msg)
{
  if ( boxesOpen )
  {
    // store message in queue (only need one element, though)
    DebugM(4,"Patch " << patchID << " proxy data arrived early, storing in buffer.\n");
    msgBuffer = msg;
    return;
  }
  DebugM(3,"Processing proxy data.\n");
  msgBuffer = NULL;
  p = msg->positionList;
  delete msg;
  positionsReady(0);
}

void ProxyPatch::receiveAll(ProxyAllMsg *msg)
{
  if ( boxesOpen )
  {
    // store message in queue (only need one element, though)
    DebugM(4,"Patch " << patchID << " proxy ALL data arrived early, storing in buffer.\n");
    msgAllBuffer = msg;
    return;
  }
  msgAllBuffer = NULL;
  DebugM(4,"Processing proxy ALL msg.\n");

  AtomMap::Object()->unregisterIDs(patchID,atomIDList);
  loadAtoms(msg->atomIDList);
  loadAtomProperties();
  AtomMap::Object()->registerIDs(patchID,msg->atomIDList);
  DebugM(4,"Processing proxy ALL msg.\n");
  p = msg->positionList;

  delete msg;

  positionsReady(1);
}

void ProxyPatch::sendResults(void)
{
  ProxyResultMsg *msg = new (MsgIndex(ProxyResultMsg)) ProxyResultMsg;
  msg->node = CMyPe();
  msg->patch = patchID;
  msg->forceList = f;
  ProxyMgr::Object()->sendResults(msg);
}

/***************************************************************************
 * RCS INFORMATION:
 *
 *	$RCSfile: ProxyPatch.C,v $
 *	$Author: ari $	$Locker:  $		$State: Exp $
 *	$Revision: 1.1001 $	$Date: 1997/02/07 05:42:32 $
 *
 ***************************************************************************
 * REVISION HISTORY:
 *
 * $Log: ProxyPatch.C,v $
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
