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
class CollectionMasterHandler;
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

#ifdef PROCTRACE_DEBUG
class DebugFileTrace;
#endif

// Instance Variables that maintain singletonness of classes

CkpvExtern(AtomMap*, AtomMap_instance);
CkpvExtern(BroadcastMgr*, BroadcastMgr_instance);
CkpvExtern(CollectionMaster*, CollectionMaster_instance);
CkpvExtern(CollectionMasterHandler*, CollectionMasterHandler_instance);
CkpvExtern(CollectionMgr*, CollectionMgr_instance);
CkpvExtern(ComputeMap*, ComputeMap_instance);
CkpvExtern(LdbCoordinator*, LdbCoordinator_instance);
CkpvExtern(Node*, Node_instance);
CkpvExtern(PatchMap*, PatchMap_instance);
CkpvExtern(PatchMgr*, PatchMgr_instance);
CkpvExtern(ProxyMgr*, ProxyMgr_instance);
CkpvExtern(ReductionMgr*, ReductionMgr_instance);
CkpvExtern(Sync*, Sync_instance);

#ifdef PROCTRACE_DEBUG
CkpvExtern(DebugFileTrace*, DebugFileTrace_instance);
#endif

// Other static Variables

CkpvExtern(PatchMgr*, PatchMap_patchMgr);
CkpvExtern(BOCgroup, BOCclass_group);
CkpvExtern(Communicate*, comm);

void ProcessorPrivateInit(void);

#endif
