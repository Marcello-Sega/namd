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


class ProxyMgr : public BOCclass
{

public:
  ProxyMgr(InitMsg *);
  ~ProxyMgr();

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
  static ProxyMgr *_instance;

};

#endif /* PATCHMGR_H */
/***************************************************************************
 * RCS INFORMATION:
 *
 *	$RCSfile: ProxyMgr.h,v $
 *	$Author: ari $	$Locker:  $		$State: Exp $
 *	$Revision: 1.1 $	$Date: 1996/12/05 21:37:53 $
 *
 ***************************************************************************
 * REVISION HISTORY:
 *
 * $Log: ProxyMgr.h,v $
 * Revision 1.1  1996/12/05 21:37:53  ari
 * Initial revision
 *
 ***************************************************************************/
