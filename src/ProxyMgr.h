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

#ifndef PROXYMGR_H
#define PROXYMGR_H

#include "ckdefs.h"
#include "chare.h"
#include "c++interface.h"

#include "main.h"
#include "NamdTypes.h"

class ProxyPatch;
class PatchMap;

class ProxyMgr : public BOCclass
{

public:
  ProxyMgr(InitMsg *);
  ~ProxyMgr();

  void removeProxies(void);
  void createProxies(void);

  void registerProxy(PatchID pid);
  void recvRegisterProxy(RegisterProxyMsg *);

  void sendResults(ProxyResultMsg *);
  void recvResults(ProxyResultMsg *);

  void sendProxyData(ProxyDataMsg *);
  void recvProxyData(ProxyDataMsg *);

  void sendProxyAtoms(ProxyAtomsMsg *);
  void recvProxyAtoms(ProxyAtomsMsg *);

  static ProxyMgr *Object() { return _instance; }
  
private:
  PatchMap *patchMap;

  static ProxyMgr *_instance;

  ProxyPatch** proxyList;
  int numProxies;

};

#endif /* PATCHMGR_H */
/***************************************************************************
 * RCS INFORMATION:
 *
 *	$RCSfile: ProxyMgr.h,v $
 *	$Author: jim $	$Locker:  $		$State: Exp $
 *	$Revision: 1.3 $	$Date: 1996/12/14 00:02:42 $
 *
 ***************************************************************************
 * REVISION HISTORY:
 *
 * $Log: ProxyMgr.h,v $
 * Revision 1.3  1996/12/14 00:02:42  jim
 * debugging ProxyAtomsMsg path to make compute creation work
 *
 * Revision 1.2  1996/12/06 03:39:09  jim
 * creation methods
 *
 * Revision 1.1  1996/12/05 21:37:53  ari
 * Initial revision
 *
 ***************************************************************************/
