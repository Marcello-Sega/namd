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


static char ident[] = "@(#)$Header: /home/cvs/namd/cvsroot/namd2/src/ProxyPatch.C,v 1.7 1996/12/14 00:02:42 jim Exp $";

#include "ckdefs.h"
#include "chare.h"
#include "c++interface.h"

#include "main.top.h"
#include "main.h"
#include "ProxyPatch.h"
#include "ProxyMgr.top.h"
#include "ProxyMgr.h"
#include "AtomMap.h"

#define MIN_DEBUG_LEVEL 3
#define  DEBUGM
#include "Debug.h"

ProxyPatch::ProxyPatch(PatchID pd) : Patch(pd)
{
  ProxyMgr::Object()->registerProxy(patchID);
}

void ProxyPatch::boxClosed(int box)
{
  if ( box == 1 )
  {
    sendResults();
  }
  if ( ! --boxesOpen )
  {
    DebugM(2,patchID << ": " << "Checking message buffer.\n");
    if ( ! msgBuffer.empty() ) receiveData(msgBuffer.pop());
  }
  else
  {
    DebugM(2,patchID << ": " << boxesOpen << " boxes left to close.\n");
  }
}

void ProxyPatch::receiveAtoms(ProxyAtomsMsg *msg)
{
  loadAtoms(*(msg->atomIDList));
  AtomMap::Object()->registerIDs(patchID,msg->atomIDList);
  delete msg;
}

void ProxyPatch::receiveData(ProxyDataMsg *msg)
{
  if ( boxesOpen )
  {
    // store message in queue
    DebugM(4,"Proxy data arrived.\n");
    msgBuffer.append(msg);
    return;
  }
  p = msg->positionList;
  delete msg;
  positionsReady();
}

void ProxyPatch::sendResults(void)
{
  ProxyResultMsg *msg = new (MsgIndex(ProxyResultMsg)) ProxyResultMsg;
  msg->patch = patchID;
  msg->forceList = f;
  ProxyMgr::Object()->sendResults(msg);
}

/***************************************************************************
 * RCS INFORMATION:
 *
 *	$RCSfile: ProxyPatch.C,v $
 *	$Author: jim $	$Locker:  $		$State: Exp $
 *	$Revision: 1.7 $	$Date: 1996/12/14 00:02:42 $
 *
 ***************************************************************************
 * REVISION HISTORY:
 *
 * $Log: ProxyPatch.C,v $
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
