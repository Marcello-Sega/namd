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

// Close all the deposit boxes
void Compute::depositAllForces() {
  PatchDepositListIter pd(patchDepositList);
  for (pd = pd.begin(); pd != pd.end(); pd++) {
    (*pd).box->close(&((*pd).f));
  }

  PatchPickupListIter pp(patchPickupList);
  for (pp = pp.begin(); pp != pp.end(); pp++) {
    (*pp).box->close(&((*pp).x));
  }

}

// Add a patch to the deposit list and register with patch
void Compute::registerForceDeposit(PatchID pid) {
  PatchDeposit pd;
  pd.pid = pid;
  pd.p = (node->patchMap).patch(pid);
  pd.box = pd.p->registerForceDeposit(cid);
  patchDepositList.load(pd);
}

// Delete a patch from the deposit list and unregister with patch
void Compute::unregisterForceDeposit(PatchID pid) {
  PatchDeposit find(pid);
  PatchDeposit* found = patchDepositList.find(find);
  found->p->unregisterForceDeposit(cid,&found->box);
  patchDepositList.del(*found);
}

// Add a patch to the pickup list and unregister with patch
void Compute::registerPositionPickup(PatchID pid) {
  PatchPickup pp;
  pp.pid = pid;
  pp.p = (node->patchMap).patch(pid);
  pp.box = pp.p->registerPositionPickup(cid);
  patchPickupList.load(pp);
  patchReadyCounter = ++numPatches;
}

// Delete a patch from the pickup list and unregister with patch
void Compute::unregisterPositionPickup(PatchID pid) {
  PatchPickup* found = patchPickupList.find(PatchPickup(pid));
  found->p->unregisterPositionPickup(cid,&found->box);
  patchPickupList.del(*found);
  patchReadyCounter = --numPatches;
}
    
Compute::~Compute() {
  // Remove the deposit registrations
  PatchDepositListIter pd(patchDepositList);
  for (pd = pd.begin(); pd != pd.end(); pd++) {
    unregisterForceDeposit((*pd).pid);
  }

  // Remove the pickup registrations
  PatchPickupListIter pp(patchPickupList);
  for (pp = pp.begin(); pp != pp.end(); pp++) {
    unregisterPositionPickup((*pp).pid);
  }
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
  CPrintf("This is the default Compute::doWork()\n");
  depositAllForces();
}

#include "Compute.bot.h"


/***************************************************************************
 * RCS INFORMATION:
 *
 *	$RCSfile: Compute.C,v $
 *	$Author: ari $	$Locker:  $		$State: Exp $
 *	$Revision: 1.1 $	$Date: 1996/10/16 08:22:39 $
 *
 ***************************************************************************
 * REVISION HISTORY:
 *
 * $Log: Compute.C,v $
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

