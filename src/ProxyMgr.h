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

class RegisterProxyMsg : public comm_object {
public:
  NodeID node;
  PatchID patch;
};

class UnregisterProxyMsg : public comm_object {
public:
  NodeID node;
  PatchID patch;
};

class ProxyAtomsMsg : public comm_object {
public:
  PatchID patch;
  AtomIDList atomIDList;
  void * pack (int *length);
  void * operator new(size_t s, int i) {return comm_object::operator new(s,i);}
  void * operator new(size_t s) { return comm_object::operator new(s); }
  void * operator new(size_t, void *ptr) { return ptr; }
  void unpack (void *in);
};

class ProxyDataMsg : public comm_object {
public:
  PatchID patch;
  PositionList positionList;
  void * pack (int *length);
  void * operator new(size_t s, int i) {return comm_object::operator new(s,i);}
  void * operator new(size_t s) { return comm_object::operator new(s); }
  void * operator new(size_t, void *ptr) { return ptr; }
  void unpack (void *in);
};

class ProxyResultMsg : public comm_object {
public:
  NodeID node;
  PatchID patch;
  ForceList forceList;
  void * pack (int *length);
  void * operator new(size_t s, int i) {return comm_object::operator new(s,i);}
  void * operator new(size_t size) { return comm_object::operator new(size); }
  void * operator new(size_t, void *ptr) { return ptr; }
  void unpack (void *in);
};

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

  void sendProxyData(ProxyDataMsg *, NodeID);
  void recvProxyData(ProxyDataMsg *);

  void sendProxyAtoms(ProxyAtomsMsg *, NodeID);
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
 *	$Revision: 1.6 $	$Date: 1996/12/17 23:58:02 $
 *
 ***************************************************************************
 * REVISION HISTORY:
 *
 * $Log: ProxyMgr.h,v $
 * Revision 1.6  1996/12/17 23:58:02  jim
 * proxy result reporting is working
 *
 * Revision 1.5  1996/12/17 17:07:41  jim
 * moved messages from main to ProxyMgr
 *
 * Revision 1.4  1996/12/17 08:56:38  jim
 * added node argument to sendProxyData and sendProxyAtoms
 *
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
