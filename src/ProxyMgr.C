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

ProxyMgr *ProxyMgr::_instance = 0;

ProxyMgr::ProxyMgr(InitMsg *) { 
  if (_instance) {
    iout << "More than one ProxyMgr!!!\n" << endi;
    Namd::die();
  }
  _instance = this;
}

ProxyMgr::~ProxyMgr() { 
  _instance = NULL;
}

void
ProxyMgr::registerProxy(PatchID pid) {
  // determine which node gets message
  NodeID node = PatchMap::Object()->node(pid);

  RegisterProxyMsg *msg = new (MsgIndex(RegisterProxyMsg)) RegisterProxyMsg;

  msg->node=CMyPe();
  msg->patch = pid;

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
  ProxyPatch *proxy = (ProxyPatch *) PatchMap::Object()->patch(msg->patch);
  proxy->receiveAtoms(msg);
}

#include "ProxyMgr.bot.h"


/***************************************************************************
 * RCS INFORMATION:
 *
 *	$RCSfile: ProxyMgr.C,v $
 *	$Author: ari $	$Locker:  $		$State: Exp $
 *	$Revision: 1.3 $	$Date: 1996/12/05 23:45:09 $
 *
 ***************************************************************************
 * REVISION HISTORY:
 *
 * $Log: ProxyMgr.C,v $
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
