/**
***  Copyright (c) 1995, 1996, 1997, 1998, 1999, 2000 by
***  The Board of Trustees of the University of Illinois.
***  All rights reserved.
**/

#include "ProcessorPrivate.h"

/*
 * Variable Definitions
 */

// Instance Variables that maintain singletonness of classes

CpvDeclare(AtomMap*, AtomMap_instance);
CpvDeclare(BroadcastMgr*, BroadcastMgr_instance);
CpvDeclare(CollectionMaster*, CollectionMaster_instance);
CpvDeclare(CollectionMgr*, CollectionMgr_instance);
CpvDeclare(ComputeMap*, ComputeMap_instance);
CpvDeclare(LJTable*, LJTable_instance);
CpvDeclare(LdbCoordinator*, LdbCoordinator_instance);
CpvDeclare(Node*, Node_instance);
CpvDeclare(PatchMap*, PatchMap_instance);
CpvDeclare(PatchMgr*, PatchMgr_instance);
CpvDeclare(ProxyMgr*, ProxyMgr_instance);
CpvDeclare(ReductionMgr*, ReductionMgr_instance);

// Other static variables

CpvDeclare(PatchMgr*, PatchMap_patchMgr);
CpvDeclare(BOCgroup, BOCclass_group);
CpvDeclare(Communicate*, comm);
CpvDeclare(Sync*, Sync_instance);

/*
 * Initialization Function to be called on every processor
 */

void ProcessorPrivateInit(void)
{
  CpvInitialize(AtomMap*, AtomMap_instance);
  CpvAccess(AtomMap_instance) = 0;
  CpvInitialize(BroadcastMgr*, BroadcastMgr_instance);
  CpvAccess(BroadcastMgr_instance) = 0;
  CpvInitialize(CollectionMaster*, CollectionMaster_instance);
  CpvAccess(CollectionMaster_instance) = 0;
  CpvInitialize(CollectionMgr*, CollectionMgr_instance);
  CpvAccess(CollectionMgr_instance) = 0;
  CpvInitialize(ComputeMap*, ComputeMap_instance);
  CpvAccess(ComputeMap_instance) = 0;
  CpvInitialize(LJTable*, LJTable_instance);
  CpvAccess(LJTable_instance) = 0;
  CpvInitialize(LdbCoordinator*, LdbCoordinator_instance);
  CpvAccess(LdbCoordinator_instance) = 0;
  CpvInitialize(Node*, Node_instance);
  CpvAccess(Node_instance) = 0;
  CpvInitialize(PatchMap*, PatchMap_instance);
  CpvAccess(PatchMap_instance) = 0;
  CpvInitialize(PatchMgr*, PatchMgr_instance);
  CpvAccess(PatchMgr_instance) = 0;
  CpvInitialize(ProxyMgr*, ProxyMgr_instance);
  CpvAccess(ProxyMgr_instance) = 0;
  CpvInitialize(ReductionMgr*, ReductionMgr_instance);
  CpvAccess(ReductionMgr_instance) = 0;
  CpvInitialize(PatchMgr*, PatchMap_patchMgr);
  CpvAccess(PatchMap_patchMgr) = 0;
  CpvInitialize(BOCgroup, BOCclass_group);
  CpvInitialize(Communicate*, comm);
  CpvAccess(comm) = 0;
  CpvInitialize(Sync*, Sync_instance);
  CpvAccess(Sync_instance) = 0;
}

