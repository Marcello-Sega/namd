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
  Patch *patch = PatchMap::Object()->patch(msg->patch);
  // patch->registerProxy(msg);
  delete msg;
}

void
ProxyMgr::sendProxyData(ProxyDataMsg *msg) {
}

void
ProxyMgr::recvProxyData(ProxyDataMsg *msg) {
}

void
ProxyMgr::sendProxyAtoms(ProxyAtomsMsg *msg) {
}

void
ProxyMgr::recvProxyAtoms(ProxyAtomsMsg *msg) {
}

#include "ProxyMgr.bot.h"


/***************************************************************************
 * RCS INFORMATION:
 *
 *	$RCSfile: ProxyMgr.C,v $
 *	$Author: ari $	$Locker:  $		$State: Exp $
 *	$Revision: 1.1 $	$Date: 1996/12/05 21:37:43 $
 *
 ***************************************************************************
 * REVISION HISTORY:
 *
 * $Log: ProxyMgr.C,v $
 * Revision 1.1  1996/12/05 21:37:43  ari
 * Initial revision
 *
 ***************************************************************************/
