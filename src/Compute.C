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
#include "WorkDistrib.h"

#include "NamdTypes.h"
#include "Templates/Box.h"
#include "Templates/OwnerBox.h"

#include "Node.h"
#include "Compute.h"

#define MIN_DEBUG_LEVEL 5
#define DEBUGM
#include "Debug.h"

void Compute::enqueueWork() {
  WorkDistrib::messageEnqueueWork(this);
}

// Signal from patch or proxy that data is ready.
// When all Patches and Proxies needed by this Compute object
// have checked-in, we are ready to enqueueWork()
void Compute::patchReady(void) { 
  if (numPatches <= 0 ) {
    DebugM(10,"Compute::patchReady()-call not valid!\n");
  } else {
    DebugM(2,"patchReadyCounter = " << patchReadyCounter << endl);
    if (! --patchReadyCounter) {
      patchReadyCounter = numPatches;
      DebugM(3,"Compute::patchReady() - enqueue()!\n");
      enqueueWork();
    }
  }
}


void Compute::doWork() {
  DebugM(5,"Default Compute::doWork() called.\n");
}

#include "Compute.bot.h"


/***************************************************************************
 * RCS INFORMATION:
 *
 *	$RCSfile: Compute.C,v $
 *	$Author: ari $	$Locker:  $		$State: Exp $
 *	$Revision: 1.777 $	$Date: 1997/01/17 19:35:32 $
 *
 ***************************************************************************
 * REVISION HISTORY:
 *
 * $Log: Compute.C,v $
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

