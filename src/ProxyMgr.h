/**
***  Copyright (c) 1995, 1996, 1997, 1998, 1999, 2000 by
***  The Board of Trustees of the University of Illinois.
***  All rights reserved.
**/

#ifndef PROXYMGR_H
#define PROXYMGR_H

#include "charm++.h"

#include "main.h"
#include "NamdTypes.h"
#include "PatchTypes.h"
#include "UniqueSet.h"
#include "UniqueSetIter.h"
#include "ProcessorPrivate.h"
#include "ProxyMgr.decl.h"

class RegisterProxyMsg : public CMessage_RegisterProxyMsg {
public:
  NodeID node;
  PatchID patch;
};

class UnregisterProxyMsg : public CMessage_UnregisterProxyMsg {
public:
  NodeID node;
  PatchID patch;
};

class ProxyAtomsMsg : public CMessage_ProxyAtomsMsg {
public:
  PatchID patch;
  AtomIDList atomIDList;
  static void* pack(ProxyAtomsMsg *msg);
  static ProxyAtomsMsg* unpack(void *ptr);
};

class ProxyDataMsg : public CMessage_ProxyDataMsg {
public:
  PatchID patch;
  Flags flags;
  PositionList positionList;
  PositionList avgPositionList;
  static void* pack(ProxyDataMsg *msg);
  static ProxyDataMsg* unpack(void *ptr);
};

class ProxyAllMsg : public CMessage_ProxyAllMsg {
public:
  PatchID patch;
  Flags flags;
  AtomIDList atomIDList;
  PositionList positionList;
  PositionList avgPositionList;
  static void* pack(ProxyAllMsg *msg);
  static ProxyAllMsg* unpack(void *ptr);
};

class ProxyResultMsg : public CMessage_ProxyResultMsg {
public:
  NodeID node;
  PatchID patch;
  ForceList forceList[Results::maxNumForces];
  static void* pack(ProxyResultMsg *msg);
  static ProxyResultMsg* unpack(void *ptr);
};

class ProxyPatch;
class PatchMap;

struct ProxyElem {
  ProxyElem() : proxyPatch(0) { };
  ProxyElem(PatchID pid) : patchID(pid), proxyPatch(0) { };
  ProxyElem(PatchID pid, ProxyPatch *p) : patchID(pid), proxyPatch(p) { };

  int hash() const { return patchID; }
  int operator==(const ProxyElem & pe) const { return patchID == pe.patchID; }

  PatchID patchID;
  ProxyPatch *proxyPatch;
};

typedef UniqueSet<ProxyElem> ProxySet;
typedef UniqueSetIter<ProxyElem> ProxySetIter;

class ProxyMgr : public BOCclass
{

public:
  ProxyMgr();
  ~ProxyMgr();

  void removeProxies(void);
  void createProxies(void);

  void createProxy(PatchID pid);
  void removeProxy(PatchID pid);

  void registerProxy(PatchID pid);
  void recvRegisterProxy(RegisterProxyMsg *);

  void sendResults(ProxyResultMsg *);
  void recvResults(ProxyResultMsg *);

  void sendProxyData(ProxyDataMsg *, NodeID);
  void recvProxyData(ProxyDataMsg *);

  void sendProxyAtoms(ProxyAtomsMsg *, NodeID);
  void recvProxyAtoms(ProxyAtomsMsg *);

  void sendProxyAll(ProxyAllMsg *, NodeID);
  void recvProxyAll(ProxyAllMsg *);

  static ProxyMgr *Object() { return CpvAccess(ProxyMgr_instance); }
  
private:

  ProxySet proxySet;
};

#endif /* PATCHMGR_H */

