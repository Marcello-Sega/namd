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

extern int proxySendSpanning, proxyRecvSpanning;

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
  CompAtomList positionList;
  CompAtomList avgPositionList;
  static void* pack(ProxyDataMsg *msg);
  static ProxyDataMsg* unpack(void *ptr);
};

class ProxyAllMsg : public CMessage_ProxyAllMsg {
public:
  PatchID patch;
  Flags flags;
  CompAtomList positionList;
  CompAtomList avgPositionList;
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

class ProxyCombinedResultMsg : public CMessage_ProxyCombinedResultMsg {
public:
  PatchID patch;
  NodeIDList nodes;
  ForceList forceList[Results::maxNumForces];
  static void* pack(ProxyCombinedResultMsg *msg);
  static ProxyCombinedResultMsg* unpack(void *ptr);
};

class ProxySpanningTreeMsg : public CMessage_ProxySpanningTreeMsg {
public:
  PatchID patch;
  NodeID  node;
  NodeIDList tree;
  static void* pack(ProxySpanningTreeMsg *msg);
  static ProxySpanningTreeMsg* unpack(void *ptr);
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
  void removeUnusedProxies(void);
  void createProxies(void);

  void createProxy(PatchID pid);
  void removeProxy(PatchID pid);

  void registerProxy(PatchID pid);
  void recvRegisterProxy(RegisterProxyMsg *);

  void unregisterProxy(PatchID pid);
  void recvUnregisterProxy(UnregisterProxyMsg *);

  void buildProxySpanningTree();
  void sendSpanningTree(ProxySpanningTreeMsg *);
  void recvSpanningTree(ProxySpanningTreeMsg *);

  void sendResults(ProxyResultMsg *);
  void recvResults(ProxyResultMsg *);
  void sendResults(ProxyCombinedResultMsg *);
  void recvResults(ProxyCombinedResultMsg *);

  void sendProxyData(ProxyDataMsg *, int, int*);
  void recvProxyData(ProxyDataMsg *);

  void sendProxyAll(ProxyAllMsg *, int, int*);
  void recvProxyAll(ProxyAllMsg *);

  static ProxyMgr *Object() { return CpvAccess(ProxyMgr_instance); }
  
  int numProxies() { return proxySet.size(); }
private:

  ProxySet proxySet;
};

#endif /* PATCHMGR_H */

