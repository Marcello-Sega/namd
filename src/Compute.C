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

#include "main.h"
#include "ckdefs.h"
#include "chare.h"
#include "c++interface.h"

#include "Compute.top.h"
#include "WorkDistrib.top.h"

#include "NamdTypes.h"
#include "Templates/Box.h"
#include "Templates/OwnerBox.h"

#include "Node.h"
#include "Compute.h"

void Compute::enqueueWork() {
  CPrintf("Compute::enqueueWork()-Sending LocalWorkMsg\n");

  LocalWorkMsg *msg = new (MsgIndex(LocalWorkMsg)) LocalWorkMsg;
  msg->compute = this; // pointer is valid since send is to local Pe
  CSendMsgBranch(WorkDistrib, enqueueWork, msg, node->workDistribGroup, 
		  CMyPe() );
}

// Signal from patch or proxy that data is ready.
// When all Patches and Proxies needed by this Compute object
// have checked-in, we are ready to enqueueWork()
void Compute::patchReady(void) { 
  if (numPatches <= 0 ) {
    CPrintf("Compute::patchReady()-call not valid!\n");
  } else {
    if (! --patchReadyCounter) {
      patchReadyCounter = numPatches;
      enqueueWork();
    }
  }
}


void Compute::doWork() {
  CPrintf("This is the default Compute::doWork() Nothing happens\n");
}

#include "Compute.bot.h"


/***************************************************************************
 * RCS INFORMATION:
 *
 *	$RCSfile: Compute.C,v $
 *	$Author: ari $	$Locker:  $		$State: Exp $
 *	$Revision: 1.2 $	$Date: 1996/10/22 19:12:16 $
 *
 ***************************************************************************
 * REVISION HISTORY:
 *
 * $Log: Compute.C,v $
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

