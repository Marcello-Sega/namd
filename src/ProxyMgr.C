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

#include "ckdefs.h"
#include "chare.h"
#include "c++interface.h"

#include "main.h"
#include "BOCgroup.h"
#include "ProxyMgr.top.h"
#include "ProxyMgr.h"
#include "Namd.h"
#include "PatchMap.inl"
#include "ProxyPatch.h"
#include "ComputeMap.h"
#include "HomePatch.h"
#include <string.h>

#include "ProcessorPrivate.h"

#define PACKDATA(LIST,TYPE,DATA) \
	memcpy((void*)DATA,(void*)(LIST.begin()),size*sizeof(TYPE))
#define UNPACKDATA(LIST,TYPE,DATA) \
	memcpy((void*)(LIST.begin()),(void*)DATA,size*sizeof(TYPE))

//#define DEBUGM
#define MIN_DEBUG_LEVEL 4
#include "Debug.h"

// Use before CSendMsg... but don't use if msg is sent to local node 
// in this case, msg will never be sent locally.  
void ProxyAtomsMsg::prepack() {
    int size = atomIDList.size();
    mylength = sizeof(PatchID) + sizeof(int) + size * sizeof(AtomID);
    char *b = mybuffer = (char*)new char[mylength];
    *((int *)b) = patch; b += sizeof(PatchID);
    *((int *)b) = size; b += sizeof(int);
    AtomID *data = (AtomID *)b;
    PACKDATA(atomIDList,AtomID,data);
}
  
  
void * ProxyAtomsMsg::pack (int *length) {
    *length = mylength;
    char *buffer = (char*)new_packbuffer(this,*length);
    memcpy(buffer, mybuffer, mylength);
    delete [] mybuffer;

    this->~ProxyAtomsMsg();
    return buffer;
}

void ProxyAtomsMsg:: unpack (void *in)
  {
    new((void*)this) ProxyAtomsMsg;
    char *buffer = (char*)in;
    patch = *((int*)buffer);
    int size = *((int*)(buffer+sizeof(int)));
    atomIDList.resize(size);
    AtomID *data = (AtomID*)(buffer+2*sizeof(int));
    UNPACKDATA(atomIDList,AtomID,data);
  }

void * ProxyDataMsg:: pack (int *length)
  {
    int size = positionList.size();
    *length = 2 * sizeof(int) + sizeof(Flags) + size * sizeof(Position);
    char *buffer = (char*)new_packbuffer(this,*length);
    *((int*)buffer) = patch;
    *((int*)(buffer+sizeof(int))) = size;
    *((Flags*)(buffer+2*sizeof(int))) = flags;
    Position *data = (Position*)(buffer+2*sizeof(int)+sizeof(Flags));
    PACKDATA(positionList,Position,data);
    this->~ProxyDataMsg();
    return buffer;
  }

void ProxyDataMsg:: unpack (void *in)
  {
    new((void*)this) ProxyDataMsg;
    char *buffer = (char*)in;
    patch = *((int*)buffer);
    int size = *((int*)(buffer+sizeof(int)));
    flags = *((Flags*)(buffer+2*sizeof(int)));
    positionList.resize(size);
    Position *data = (Position*)(buffer+2*sizeof(int)+sizeof(Flags));
    UNPACKDATA(positionList,Position,data);
  }

void * ProxyAllMsg:: pack (int *length)
  {
    int size = positionList.size();
    if (size != atomIDList.size()) {
      iout << "ProxyAllMsg::pack() - Bad News, sizes don't match!" << endi;
    }


    *length = 2 * sizeof(int) + sizeof(Flags) + size * sizeof(Position) + size * sizeof(AtomID);
    char *buffer = (char*)new_packbuffer(this,*length);
    *((int*)buffer) = patch;
    *((int*)(buffer+sizeof(int))) = size;
    *((Flags*)(buffer+2*sizeof(int))) = flags;
    Position *data = (Position*)(buffer+2*sizeof(int)+sizeof(Flags));
    PACKDATA(positionList,Position,data);
    AtomID *data2 = (AtomID*)(buffer+2*sizeof(int)+sizeof(Flags)+size*sizeof(Position));
    PACKDATA(atomIDList,AtomID,data2);
    this->~ProxyAllMsg();
    return buffer;
  }

void ProxyAllMsg:: unpack (void *in)
  {
    new((void*)this) ProxyAllMsg;
    char *buffer = (char*)in;
    patch = *((int*)buffer);
    int size = *((int*)(buffer+sizeof(int)));
    flags = *((Flags*)(buffer+2*sizeof(int)));
    positionList.resize(size);
    Position *data = (Position*)(buffer+2*sizeof(int)+sizeof(Flags));
    UNPACKDATA(positionList,Position,data);
    atomIDList.resize(size);
    AtomID *data2 = (AtomID*)(buffer+2*sizeof(int)+sizeof(Flags)+size*sizeof(Position));
    UNPACKDATA(atomIDList,AtomID,data2);
  }

void * ProxyResultMsg:: pack (int *length)
  {
    *length = ( 4 + Results::maxNumForces ) * sizeof(int);
    for ( int j = 0; j < Results::maxNumForces; ++j )
    {
      *length += sizeof(Force) * forceList[j].size();
    }
    char *buffer = (char*)new_packbuffer(this,*length);
    char *b = buffer;
    *((int*)b) = node;  b += sizeof(int);
    *((int*)b) = patch;  b += sizeof(int);
    for ( j = 0; j < Results::maxNumForces; ++j )
    {
      int size = forceList[j].size();
      memcpy((void*)b,(void*)(&size),sizeof(int));  b += sizeof(int);
      memcpy((void*)b,(void*)(forceList[j].begin()),size*sizeof(Force));
      b += size*sizeof(Force);
    }
    this->~ProxyResultMsg();
    return buffer;
  }

void ProxyResultMsg:: unpack (void *in)
  {
    new((void*)this) ProxyResultMsg;
    char *b = (char*)in;
    node = *((int*)b);  b += sizeof(int);
    patch = *((int*)b);  b += sizeof(int);
    for ( int j = 0; j < Results::maxNumForces; ++j )
    {
      int size;
      memcpy((void*)(&size),(void*)b,sizeof(int));  b += sizeof(int);
      forceList[j].resize(size);
      memcpy((void*)(forceList[j].begin()),(void*)b,size*sizeof(Force));
      b += size*sizeof(Force);
    }
  }

ProxyMgr::ProxyMgr(InitMsg *) { 
  if (CpvAccess(ProxyMgr_instance)) {
    Namd::die();
  }
  CpvAccess(ProxyMgr_instance) = this;
}

ProxyMgr::~ProxyMgr() { 
  removeProxies();
  CpvAccess(ProxyMgr_instance) = NULL;
}

void ProxyMgr::removeProxies(void)
{
  ProxySetIter pi(proxySet);
  for ( pi = pi.begin(); pi != pi.end(); pi++)
  {
    delete pi->proxyPatch;
  }
  proxySet.clear();
}

// Figure out which proxies we need and create them
void ProxyMgr::createProxies(void)
{
  // Delete the old proxies.
  removeProxies();

  PatchMap *patchMap = PatchMap::Object();
  int numPatches = patchMap->numPatches();
  int myNode = CMyPe();
  enum PatchFlag { Unknown, Home, NeedProxy };
  int *patchFlag = new int[numPatches]; 
  int i, j;

  // Note all home patches.
  for ( i = 0; i < numPatches; ++i )
  {
    patchFlag[i] = ( patchMap->node(i) == myNode ) ? Home : Unknown;
  }

  // Add all upstream neighbors.
  PatchID neighbors[PatchMap::MaxOneAway];
  for ( i = 0; i < numPatches; ++i )
  {
    if ( patchMap->node(i) != myNode ) 
      continue;
    int numNeighbors = patchMap->upstreamNeighbors(i,neighbors);
    for ( j = 0; j < numNeighbors; ++j )
    {
      if ( ! patchFlag[neighbors[j]] ) {
	patchFlag[neighbors[j]] = NeedProxy;
      }
    }
  }

  // Check all patch-based compute objects.
  ComputeMap *computeMap = ComputeMap::Object();
  int nc = computeMap->numComputes();
  for ( i = 0; i < nc; ++i )
  {
    if ( computeMap->node(i) != myNode || !computeMap->isPatchBased(i) ) 
      continue;
    int numPid = computeMap->numPids(i);
    for ( j = 0; j < numPid; ++j )
    {
      int pid = computeMap->pid(i,j);
      if ( ! patchFlag[pid] ) {
	patchFlag[pid] = NeedProxy;
      }
    }
  }
  
  // Create proxy list
  for ( i = 0; i < numPatches; ++i ) {
    if ( patchFlag[i] == NeedProxy )
    { // create proxy patch
      ProxyPatch *proxy = new ProxyPatch(i);
      proxySet.add(ProxyElem(i, proxy));
      patchMap->registerPatch(i, proxy);
    }
  }
  delete[] patchFlag;
}

void
ProxyMgr::createProxy(PatchID pid) {
  Patch *p = PatchMap::Object()->patch(pid);
  if (!p) {
     DebugM(4,"createProxy("<<pid<<")\n");
     ProxyPatch *proxy = new ProxyPatch(pid);
     proxySet.add(ProxyElem(pid,proxy));
     PatchMap::Object()->registerPatch(pid,proxy);
  }
  else {
     DebugM(4,"createProxy("<<pid<<") found " << p->getPatchID() << "\n");
  }
    
}

void
ProxyMgr::removeProxy(PatchID pid) {
  ProxyElem *p = proxySet.find(ProxyElem(pid));
  if (p) { 
    delete p->proxyPatch;
    proxySet.del(ProxyElem(pid));
  }
}
  
void
ProxyMgr::registerProxy(PatchID pid) {
  // determine which node gets message
  NodeID node = PatchMap::Object()->node(pid);

  RegisterProxyMsg *msg = new (MsgIndex(RegisterProxyMsg)) RegisterProxyMsg;

  msg->node=CMyPe();
  msg->patch = pid;

  CSendMsgBranch(ProxyMgr, recvRegisterProxy, RegisterProxyMsg, msg, CpvAccess(BOCclass_group).proxyMgr, node);
}

void
ProxyMgr::recvRegisterProxy(RegisterProxyMsg *msg) {
  HomePatch *homePatch = (HomePatch *)PatchMap::Object()->patch(msg->patch);
  homePatch->registerProxy(msg); // message deleted in registerProxy()
}

void
ProxyMgr::sendResults(ProxyResultMsg *msg) {
  NodeID node = PatchMap::Object()->node(msg->patch);
  CSendMsgBranch(ProxyMgr, recvResults, ProxyResultMsg, msg, CpvAccess(BOCclass_group).proxyMgr, node);
}

void
ProxyMgr::recvResults(ProxyResultMsg *msg) {
  HomePatch *home = (HomePatch *) PatchMap::Object()->patch(msg->patch);
  home->receiveResults(msg); // delete done in HomePatch::receiveResults()
}

void
ProxyMgr::sendProxyData(ProxyDataMsg *msg, NodeID node) {
  CSendMsgBranch(ProxyMgr, recvProxyData, ProxyDataMsg, msg, CpvAccess(BOCclass_group).proxyMgr, node);
}

void
ProxyMgr::recvProxyData(ProxyDataMsg *msg) {
  ProxyPatch *proxy = (ProxyPatch *) PatchMap::Object()->patch(msg->patch);
  proxy->receiveData(msg); // deleted in ProxyPatch::receiveAtoms()
}

void
ProxyMgr::sendProxyAtoms(ProxyAtomsMsg *msg, NodeID node) {
  CSendMsgBranch(ProxyMgr, recvProxyAtoms, ProxyAtomsMsg, msg, CpvAccess(BOCclass_group).proxyMgr, node);
}

void
ProxyMgr::recvProxyAtoms(ProxyAtomsMsg *msg) {
  ProxyPatch *proxy = (ProxyPatch *) PatchMap::Object()->patch(msg->patch);
  proxy->receiveAtoms(msg); // deleted in ProxyPatch::receiveAtoms()
}

void
ProxyMgr::sendProxyAll(ProxyAllMsg *msg, NodeID node) {
  CSendMsgBranch(ProxyMgr, recvProxyAll, ProxyAllMsg, msg, CpvAccess(BOCclass_group).proxyMgr, node);
}

void
ProxyMgr::recvProxyAll(ProxyAllMsg *msg) {
  ProxyPatch *proxy = (ProxyPatch *) PatchMap::Object()->patch(msg->patch);
  proxy->receiveAll(msg); // delete done in ProxyPatch::receiveAll()
}

#include "ProxyMgr.bot.h"


/***************************************************************************
 * RCS INFORMATION:
 *
 *	$RCSfile: ProxyMgr.C,v $
 *	$Author: milind $	$Locker:  $		$State: Exp $
 *	$Revision: 1.1021 $	$Date: 1998/02/10 23:30:31 $
 *
 ***************************************************************************
 * REVISION HISTORY:
 *
 * $Log: ProxyMgr.C,v $
 * Revision 1.1021  1998/02/10 23:30:31  milind
 * Fixed to reflect the current changes to Charm++ translator.
 *
 * Revision 1.1020  1997/12/22 21:29:26  jim
 * Proxies no longer send empty arrays back to HomePatch.  Requires some new
 * flags to be set correctly in Sequencer in order to work.  These are:
 *   maxForceMerged - this and faster are added into Results::normal array
 *   maxForceUsed - all forces slower than this are discarded (assumed zero)
 * Generally maxForceMerged doesn't change but maxForceUsed depends on timestep.
 *
 * Revision 1.1019  1997/12/19 23:42:36  jim
 * Replaced assignments with memcpys and reordered memcpys for efficiency.
 *
 * Revision 1.1018  1997/12/10 17:53:35  milind
 * Removed the dcd file already exists error. Now, if a dcd file already exists,
 * it is moved to a .bak before writing new dcd file.
 *
 * Revision 1.1017  1997/11/07 20:17:46  milind
 * Made NAMD to run on shared memory machines.
 *
 * Revision 1.1016  1997/10/06 00:12:35  jim
 * Added PatchMap.inl, sped up cycle-boundary tuple code.
 *
 * Revision 1.1015  1997/09/28 22:36:53  jim
 * Modified tuple-based computations to not duplicate calculations and
 * only require "upstream" proxies.
 *
 * Revision 1.1014  1997/04/10 09:14:09  ari
 * Final debugging for compute migration / proxy creation for load balancing.
 * Lots of debug code added, mostly turned off now.
 * Fixed bug in PositionBox when Patch had no dependencies.
 * Eliminated use of cout and misuse of iout in numerous places.
 *                                            Ari & Jim
 *
 * Revision 1.1013  1997/04/08 07:08:55  ari
 * Modification for dynamic loadbalancing - moving computes
 * Still bug in new computes or usage of proxies/homepatches.
 * Works if ldbStrategy is none as before.
 *
 * Revision 1.1012  1997/03/20 23:53:48  ari
 * Some changes for comments. Copyright date additions.
 * Hooks for base level update of Compute objects from ComputeMap
 * by ComputeMgr.  Useful for new compute migration functionality.
 *
 * Revision 1.1011  1997/03/12 22:06:46  jim
 * First step towards multiple force returns and multiple time stepping.
 *
 * Revision 1.1010  1997/02/28 04:47:11  jim
 * Full electrostatics now works with fulldirect on one node.
 *
 * Revision 1.1009  1997/02/26 16:53:16  ari
 * Cleaning and debuging for memory leaks.
 * Adding comments.
 * Removed some dead code due to use of Quiescense detection.
 *
 * Revision 1.1008  1997/02/25 18:26:22  nealk
 * More DPMTA debugging
 * Disabled ProxyMgr debugging
 *
 * Revision 1.1007  1997/02/17 23:47:05  ari
 * Added files for cleaning up atom migration code
 *
 * Revision 1.1006  1997/02/13 23:17:19  ari
 * Fixed a final bug in AtomMigration - numatoms in ComputePatchPair.C not
 * set correctly in atomUpdate()
 *
 * Revision 1.1005  1997/02/13 16:17:18  ari
 * Intermediate debuging commit - working to fix deep bug in migration?
 *
 * Revision 1.1004  1997/02/13 04:43:14  jim
 * Fixed initial hanging (bug in PatchMap, but it still shouldn't have
 * happened) and saved migration messages in the buffer from being
 * deleted, but migration still dies (even on one node).
 *
 * Revision 1.1003  1997/02/11 18:51:55  ari
 * Modified with #ifdef DPMTA to safely eliminate DPMTA codes
 * fixed non-buffering of migration msgs
 * Migration works on multiple processors
 *
 * Revision 1.1002  1997/02/07 17:39:40  ari
 * More debugging for atomMigration.
 * Using -w on CC got us some minor fixes
 * using purify got us a major memory problem due to bad sizing of dummy force
 *
 * Revision 1.1001  1997/02/06 21:56:34  jim
 * Fixed bugs.
 *
 * Revision 1.1000  1997/02/06 15:59:11  ari
 * Resetting CVS to merge branches back into the main trunk.
 * We will stick to main trunk development as suggested by CVS manual.
 * We will set up tags to track fixed points of development/release
 * as suggested by CVS manual - all praise the CVS manual.
 *
 * Revision 1.779  1997/02/06 15:53:26  ari
 * Updating Revision Line, getting rid of branches
 *
 * Revision 1.778.2.1  1997/02/05 22:18:21  ari
 * Added migration code - Currently the framework is
 * there with compiling code.  This version does
 * crash shortly after migration is complete.
 * Migration appears to complete, but Patches do
 * not appear to be left in a correct state.
 *
 * Revision 1.778  1997/01/28 00:31:17  ari
 * internal release uplevel to 1.778
 *
 * Revision 1.777.2.1  1997/01/27 22:45:36  ari
 * Basic Atom Migration Code added.
 * Added correct magic first line to .h files for xemacs to go to C++ mode.
 * Compiles and runs without migration turned on.
 *
 * Revision 1.777  1997/01/17 19:36:52  ari
 * Internal CVS leveling release.  Start development code work
 * at 1.777.1.1.
 *
 * Revision 1.12  1997/01/15 17:09:43  ari
 * minor changes
 *
 * Revision 1.11  1996/12/17 23:58:02  jim
 * proxy result reporting is working
 *
 * Revision 1.10  1996/12/17 17:07:41  jim
 * moved messages from main to ProxyMgr
 *
 * Revision 1.9  1996/12/17 08:56:38  jim
 * added node argument to sendProxyData and sendProxyAtoms
 *
 * Revision 1.8  1996/12/14 00:02:42  jim
 * debugging ProxyAtomsMsg path to make compute creation work
 *
 * Revision 1.7  1996/12/13 19:39:55  jim
 * added debugging, looking for error in PatchMap sending
 *
 * Revision 1.6  1996/12/13 08:54:53  jim
 * fixed initialization bug
 *
 * Revision 1.5  1996/12/10 00:13:12  ari
 * *** empty log message ***
 *
 * Revision 1.4  1996/12/06 03:39:09  jim
 * creation methods
 *
 * Revision 1.3  1996/12/05 23:45:09  ari
 * *** empty log message ***
 *
 * Revision 1.2  1996/12/05 22:17:38  jim
 * added functions
 *
 * Revision 1.1  1996/12/05 21:37:43  ari
 * Initial revision
 *
 ***************************************************************************/
