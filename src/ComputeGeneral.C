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

#include "NamdTypes.h"
#include "Templates/Box.h"
#include "Templates/OwnerBox.h"

#include "Node.h"
#include "Patch.h"
#include "ComputeGeneral.h"


// Close all the deposit boxes
void ComputeGeneral::depositAllForces() {
  PatchDepositListIter pd(patchDepositList);
  for (pd = pd.begin(); pd != pd.end(); pd++) {
    (*pd).box->close(&((*pd).r));
  }

  PatchPickupListIter pp(patchPickupList);
  for (pp = pp.begin(); pp != pp.end(); pp++) {
    (*pp).box->close(&((*pp).x));
  }

}

// Add a patch to the deposit list and register with patch
void ComputeGeneral::registerForceDeposit(PatchID pid) {
  PatchDeposit pd;
  pd.pid = pid;
  pd.p = patchMap->patch(pid);
  pd.box = pd.p->registerForceDeposit(cid);
  patchDepositList.load(pd);
}

// Delete a patch from the deposit list and unregister with patch
void ComputeGeneral::unregisterForceDeposit(PatchID pid) {
  PatchDeposit find(pid);
  PatchDeposit* found = patchDepositList.find(find);
  found->p->unregisterForceDeposit(cid,&found->box);
  patchDepositList.del(*found);
}

// Add a patch to the pickup list and unregister with patch
void ComputeGeneral::registerPositionPickup(PatchID pid) {
  PatchPickup pp;
  pp.pid = pid;
  pp.p = patchMap->patch(pid);
  pp.box = pp.p->registerPositionPickup(cid);
  patchPickupList.load(pp);
  setNumPatches(getNumPatches() + 1);
}

// Delete a patch from the pickup list and unregister with patch
void ComputeGeneral::unregisterPositionPickup(PatchID pid) {
  PatchPickup* found = patchPickupList.find(PatchPickup(pid));
  found->p->unregisterPositionPickup(cid,&found->box);
  patchPickupList.del(*found);
  setNumPatches(getNumPatches() - 1);
}
    
ComputeGeneral::~ComputeGeneral() {
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


void ComputeGeneral::doWork() {
  CPrintf("This is the default ComputeGeneral::doWork()\n");
  depositAllForces();
}


/***************************************************************************
 * RCS INFORMATION:
 *
 *	$RCSfile: ComputeGeneral.C,v $
 *	$Author: jim $	$Locker:  $		$State: Exp $
 *	$Revision: 1.1002 $	$Date: 1997/03/18 18:08:51 $
 *
 ***************************************************************************
 * REVISION HISTORY:
 *
 * $Log: ComputeGeneral.C,v $
 * Revision 1.1002  1997/03/18 18:08:51  jim
 * Revamped collection system to ensure ordering and eliminate
 * unnecessary collections.  Also reduced make dependencies.
 *
 * Revision 1.1001  1997/03/12 22:06:33  jim
 * First step towards multiple force returns and multiple time stepping.
 *
 * Revision 1.1000  1997/02/06 15:57:55  ari
 * Resetting CVS to merge branches back into the main trunk.
 * We will stick to main trunk development as suggested by CVS manual.
 * We will set up tags to track fixed points of development/release
 * as suggested by CVS manual - all praise the CVS manual.
 *
 * Revision 1.778  1997/01/28 00:30:08  ari
 * internal release uplevel to 1.778
 *
 * Revision 1.777  1997/01/17 19:35:41  ari
 * Internal CVS leveling release.  Start development code work
 * at 1.777.1.1.
 *
 * Revision 1.2  1996/10/29 23:35:27  ari
 * *** empty log message ***
 *
 * Revision 1.1  1996/10/22 19:15:14  ari
 * Initial revision
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

