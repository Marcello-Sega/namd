/**
***  Copyright (c) 1995, 1996, 1997, 1998, 1999, 2000 by
***  The Board of Trustees of the University of Illinois.
***  All rights reserved.
**/

#include "charm++.h"

#ifndef PROCESSOR_PRIVATE_H
#define PROCESSOR_PRIVATE_H

#include "BOCgroup.h"

class AtomMap;
class BroadcastMgr;
class CollectionMaster;
class CollectionMgr;
class ComputeMap;
class LdbCoordinator;
class Node;
class PatchMap;
class PatchMgr;
class ProxyMgr;
class ReductionMgr;
class Communicate;
class Sync;

// Instance Variables that maintain singletonness of classes

CpvExtern(AtomMap*, AtomMap_instance);
CpvExtern(BroadcastMgr*, BroadcastMgr_instance);
CpvExtern(CollectionMaster*, CollectionMaster_instance);
CpvExtern(CollectionMgr*, CollectionMgr_instance);
CpvExtern(ComputeMap*, ComputeMap_instance);
CpvExtern(LdbCoordinator*, LdbCoordinator_instance);
CpvExtern(Node*, Node_instance);
CpvExtern(PatchMap*, PatchMap_instance);
CpvExtern(PatchMgr*, PatchMgr_instance);
CpvExtern(ProxyMgr*, ProxyMgr_instance);
CpvExtern(ReductionMgr*, ReductionMgr_instance);
CpvExtern(Sync*, Sync_instance);

// Other static Variables

CpvExtern(PatchMgr*, PatchMap_patchMgr);
CpvExtern(BOCgroup, BOCclass_group);
CpvExtern(Communicate*, comm);

void ProcessorPrivateInit(void);

#endif
