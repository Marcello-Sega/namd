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
#include "PatchMap.h"
#include "ProxyPatch.h"
#include "ComputeMap.h"
#include "HomePatch.h"

#define DEBUGM
#define MIN_DEBUG_LEVEL 4
#include "Debug.h"

void * ProxyAtomsMsg::pack (int *length)
  {
    int size = atomIDList.size();
    DebugM(3,"Size of atomIDList " << size << "\n");
    *length = 2 * sizeof(int) + size * sizeof(AtomID);
    char *buffer = (char*)new_packbuffer(this,*length);
    *((int*)buffer) = patch;
    *((int*)(buffer+sizeof(int))) = size;
    AtomID *data = (AtomID*)(buffer+2*sizeof(int));
    for ( int i = 0; i < size; ++i )
      data[i] = atomIDList[i];
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
    for ( int i = 0; i < size; ++i )
      atomIDList[i] = data[i];
  }

void * ProxyDataMsg:: pack (int *length)
  {
    int size = positionList.size();
    *length = 2 * sizeof(int) + size * sizeof(Position);
    char *buffer = (char*)new_packbuffer(this,*length);
    *((int*)buffer) = patch;
    *((int*)(buffer+sizeof(int))) = size;
    Position *data = (Position*)(buffer+2*sizeof(int));
    for ( int i = 0; i < size; ++i )
      data[i] = positionList[i];
    this->~ProxyDataMsg();
    return buffer;
  }

void ProxyDataMsg:: unpack (void *in)
  {
    new((void*)this) ProxyDataMsg;
    char *buffer = (char*)in;
    patch = *((int*)buffer);
    int size = *((int*)(buffer+sizeof(int)));
    positionList.resize(size);
    Position *data = (Position*)(buffer+2*sizeof(int));
    for ( int i = 0; i < size; ++i )
      positionList[i] = data[i];
  }

void * ProxyAllMsg:: pack (int *length)
  {
    int i;

    int size = positionList.size();
    if (size != atomIDList.size()) {
      iout << "ProxyAllMsg::pack() - Bad News, sizes don't match!" << endi;
    }


    *length = 2 * sizeof(int) + size * sizeof(Position) + size * sizeof(AtomID);
    char *buffer = (char*)new_packbuffer(this,*length);
    *((int*)buffer) = patch;
    *((int*)(buffer+sizeof(int))) = size;
    Position *data = (Position*)(buffer+2*sizeof(int));
    for ( i = 0; i < size; ++i )
      data[i] = positionList[i];
    AtomID *data2 = (AtomID*)(buffer+2*sizeof(int)+size*sizeof(Position));
    for ( i = 0; i < size; ++i )
      data2[i] = atomIDList[i];
    this->~ProxyAllMsg();
    return buffer;
  }

void ProxyAllMsg:: unpack (void *in)
  {
    int i;

    new((void*)this) ProxyAllMsg;
    char *buffer = (char*)in;
    patch = *((int*)buffer);
    int size = *((int*)(buffer+sizeof(int)));
    positionList.resize(size);
    Position *data = (Position*)(buffer+2*sizeof(int));
    for ( i = 0; i < size; ++i )
      positionList[i] = data[i];
    atomIDList.resize(size);
    AtomID *data2 = (AtomID*)(buffer+2*sizeof(int)+size*sizeof(Position));
    for ( i = 0; i < size; ++i )
      atomIDList[i] = data2[i];
  }

void * ProxyResultMsg:: pack (int *length)
  {
    int size = forceList.size();
    *length = 4 * sizeof(int) + size * sizeof(Force);
    char *buffer = (char*)new_packbuffer(this,*length);
    *((int*)buffer) = node;
    *((int*)(buffer+sizeof(int))) = patch;
    *((int*)(buffer+2*sizeof(int))) = size;
    Force *data = (Force*)(buffer+4*sizeof(int));
    for ( int i = 0; i < size; ++i )
      data[i] = forceList[i];
    this->~ProxyResultMsg();
    return buffer;
  }

void ProxyResultMsg:: unpack (void *in)
  {
    new((void*)this) ProxyResultMsg;
    char *buffer = (char*)in;
    node = *((int*)buffer);
    patch = *((int*)(buffer+sizeof(int)));
    int size = *((int*)(buffer+2*sizeof(int)));
    forceList.resize(size);
    Force *data = (Force*)(buffer+4*sizeof(int));
    for ( int i = 0; i < size; ++i )
      forceList[i] = data[i];
  }

ProxyMgr *ProxyMgr::_instance = 0;

ProxyMgr::ProxyMgr(InitMsg *) : numProxies(0), proxyList(NULL) { 
  DebugM(1, "::ProxyMgr() - my pe is " << CMyPe() << endl );
  if (_instance) {
    iout << "More than one ProxyMgr!!!\n" << endi;
    // Namd::die();
  }
  _instance = this;

  patchMap = PatchMap::Instance();
}

ProxyMgr::~ProxyMgr() { 
  _instance = NULL;
}

void ProxyMgr::removeProxies(void)
{
  for ( int i = 0; i < numProxies; ++i )
  {
    delete proxyList[i];
  }
  delete [] proxyList;
  proxyList = NULL;
  numProxies = 0;
}

void ProxyMgr::createProxies(void)
{
  // Delete the old proxies.
  removeProxies();

  // Figure out which proxies we will be needing.
  PatchMap *pmap = PatchMap::Object();
  int n = pmap->numPatches();
  int myNode = CMyPe();
  int *pflags = new int[n]; // 0 = unknown, 1 = home, 2 = proxy needed
  int i, j;
  // Note all home patches.
  for ( i = 0; i < n; ++i )
  {
    pflags[i] = pmap->node(i) == myNode ? 1 : 0;
  }
  // Check all two-away neighbors.
  PatchID neighbors[PatchMap::MaxOneAway + PatchMap::MaxTwoAway];
  for ( i = 0; i < n; ++i )
  {
    int nn = pmap->oneAwayNeighbors(i,neighbors);
    nn += pmap->twoAwayNeighbors(i,neighbors+nn);
    for ( j = 0; j < nn; ++j )
    {
      if ( ! pflags[neighbors[j]] ) pflags[neighbors[j]] = 2;
    }
  }
  // Check all patch-based compute objects.
  ComputeMap *cmap = ComputeMap::Object();
  int nc = cmap->numComputes();
  for ( i = 0; i < nc; ++i )
  {
    if ( cmap->node(i) != myNode || ! cmap->isPatchBased(i) ) continue;
    int ncp = cmap->numPids(i);
    for ( j = 0; j < ncp; ++j )
    {
      int pid = cmap->pid(i,j);
      if ( ! pflags[pid] ) pflags[pid] = 2;
    }
  }
  
  // Create proxies.
  for ( i = 0; i < n; ++i ) if ( pflags[i] == 2 ) ++numProxies;
  proxyList = new ProxyPatch*[numProxies];
  numProxies = 0;
  for ( i = 0; i < n; ++i ) if ( pflags[i] == 2 )
  {
    ProxyPatch *proxy = new ProxyPatch(i);
    proxyList[numProxies++] = proxy;
    patchMap->registerPatch(i, proxy);
  }
}

void
ProxyMgr::registerProxy(PatchID pid) {
  // determine which node gets message
  NodeID node = PatchMap::Object()->node(pid);

  RegisterProxyMsg *msg = new (MsgIndex(RegisterProxyMsg)) RegisterProxyMsg;

  msg->node=CMyPe();
  msg->patch = pid;

  DebugM(1,"For patch " << pid << " registering proxy on node " << CMyPe() << " with home patch on node " << node << endl);

  CSendMsgBranch(ProxyMgr, recvRegisterProxy, msg, group.proxyMgr, node);
}

void
ProxyMgr::recvRegisterProxy(RegisterProxyMsg *msg) {
  DebugM(1,"For patch " << msg->patch << " registering proxy on node " << msg->node << " with home patch on node " << CMyPe() << endl);
  HomePatch *homePatch = (HomePatch *)PatchMap::Object()->patch(msg->patch);
  homePatch->registerProxy(msg);
  delete msg;
}

void
ProxyMgr::sendResults(ProxyResultMsg *msg) {
  NodeID node = PatchMap::Object()->node(msg->patch);
  DebugM(1,"For patch " << msg->patch << " sending results from proxy on node " << CMyPe() << " to home patch on node " << node << endl);
  CSendMsgBranch(ProxyMgr, recvResults, msg, group.proxyMgr, node);
}

void
ProxyMgr::recvResults(ProxyResultMsg *msg) {
  DebugM(1,"For patch " << msg->patch << " received results from proxy on node " << msg->node << endl);
  HomePatch *home = (HomePatch *) PatchMap::Object()->patch(msg->patch);
  home->receiveResults(msg);
  delete msg;
}

void
ProxyMgr::sendProxyData(ProxyDataMsg *msg, NodeID node) {
  DebugM(1,"For patch " << msg->patch << " sending data for proxy on node " << node << " from home patch on node " << CMyPe() << endl);
  CSendMsgBranch(ProxyMgr, recvProxyData, msg, group.proxyMgr, node);
}

void
ProxyMgr::recvProxyData(ProxyDataMsg *msg) {
  DebugM(1,"For patch " << msg->patch << " received data for proxy on node " << CMyPe() << endl);
  ProxyPatch *proxy = (ProxyPatch *) PatchMap::Object()->patch(msg->patch);
  proxy->receiveData(msg);
}

void
ProxyMgr::sendProxyAtoms(ProxyAtomsMsg *msg, NodeID node) {
  DebugM(1,"For patch " << msg->patch << " sending atoms for proxy on node " << node << " from home patch on node " << CMyPe() << endl);
  CSendMsgBranch(ProxyMgr, recvProxyAtoms, msg, group.proxyMgr, node);
}

void
ProxyMgr::recvProxyAtoms(ProxyAtomsMsg *msg) {
  DebugM(1,"For patch " << msg->patch << " received atoms for proxy on node " << CMyPe() << endl);
  ProxyPatch *proxy = (ProxyPatch *) PatchMap::Object()->patch(msg->patch);
  proxy->receiveAtoms(msg);
}

void
ProxyMgr::sendProxyAll(ProxyAllMsg *msg, NodeID node) {
  DebugM(1,"For patch " << msg->patch << " sending data for proxy on node " << node << " from home patch on node " << CMyPe() << endl);
  CSendMsgBranch(ProxyMgr, recvProxyAll, msg, group.proxyMgr, node);
}

void
ProxyMgr::recvProxyAll(ProxyAllMsg *msg) {
  DebugM(1,"For patch " << msg->patch << " received data for proxy on node " << CMyPe() << endl);
  ProxyPatch *proxy = (ProxyPatch *) PatchMap::Object()->patch(msg->patch);
  proxy->receiveAll(msg);
}

#include "ProxyMgr.bot.h"


/***************************************************************************
 * RCS INFORMATION:
 *
 *	$RCSfile: ProxyMgr.C,v $
 *	$Author: ari $	$Locker:  $		$State: Exp $
 *	$Revision: 1.1002 $	$Date: 1997/02/07 17:39:40 $
 *
 ***************************************************************************
 * REVISION HISTORY:
 *
 * $Log: ProxyMgr.C,v $
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
