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
#include "charm++.h"

#include "NamdTypes.h"
#include "Box.h"
#include "OwnerBox.h"

#include "Node.h"
#include "Patch.h"
#include "PatchMap.h"
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
  CkPrintf("This is the default ComputeGeneral::doWork()\n");
  depositAllForces();
}

