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

#define DEBUGM
#include "Debug.h"

ProxyMgr *ProxyMgr::_instance = 0;

ProxyMgr::ProxyMgr(InitMsg *) : numProxies(0), proxyList(NULL) { 
  DebugM(1, "::ProxyMgr() - my pe is " << CMyPe() << endl );
  if (_instance) {
    iout << "More than one ProxyMgr!!!\n" << endi;
    // Namd::die();
  }
  _instance = this;
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
    proxyList[numProxies++] = new ProxyPatch(i);
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
  HomePatch *homePatch = (HomePatch *)PatchMap::Object()->patch(msg->patch);
  homePatch->registerProxy(msg);
  delete msg;
}

void
ProxyMgr::sendResults(ProxyResultMsg *msg) {
}

void
ProxyMgr::recvResults(ProxyResultMsg *msg) {
}

void
ProxyMgr::sendProxyData(ProxyDataMsg *msg) {
  NodeID node = PatchMap::Object()->node(msg->patch);
  CSendMsgBranch(ProxyMgr, recvProxyData, msg, group.proxyMgr, node);
}

void
ProxyMgr::recvProxyData(ProxyDataMsg *msg) {
  DebugM(1,"For patch " << msg->patch << " received data for proxy on node " << CMyPe() << endl);
  ProxyPatch *proxy = (ProxyPatch *) PatchMap::Object()->patch(msg->patch);
  proxy->receiveData(msg);
}

void
ProxyMgr::sendProxyAtoms(ProxyAtomsMsg *msg) {
  NodeID node = PatchMap::Object()->node(msg->patch);
  CSendMsgBranch(ProxyMgr, recvProxyAtoms, msg, group.proxyMgr, node);
}

void
ProxyMgr::recvProxyAtoms(ProxyAtomsMsg *msg) {
  DebugM(1,"For patch " << msg->patch << " received atoms for proxy on node " << CMyPe() << endl);
  ProxyPatch *proxy = (ProxyPatch *) PatchMap::Object()->patch(msg->patch);
  proxy->receiveAtoms(msg);
}

#include "ProxyMgr.bot.h"


/***************************************************************************
 * RCS INFORMATION:
 *
 *	$RCSfile: ProxyMgr.C,v $
 *	$Author: jim $	$Locker:  $		$State: Exp $
 *	$Revision: 1.7 $	$Date: 1996/12/13 19:39:55 $
 *
 ***************************************************************************
 * REVISION HISTORY:
 *
 * $Log: ProxyMgr.C,v $
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
