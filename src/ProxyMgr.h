//-*-c++-*-
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

class ProxyAllMsg : public comm_object {
public:
  PatchID patch;
  AtomIDList atomIDList;
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

  void sendProxyAll(ProxyAllMsg *, NodeID);
  void recvProxyAll(ProxyAllMsg *);

  static ProxyMgr *Object() { return _instance; }
  
private:
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
 *	$Revision: 1.1001 $	$Date: 1997/02/13 04:43:15 $
 *
 ***************************************************************************
 * REVISION HISTORY:
 *
 * $Log: ProxyMgr.h,v $
 * Revision 1.1001  1997/02/13 04:43:15  jim
 * Fixed initial hanging (bug in PatchMap, but it still shouldn't have
 * happened) and saved migration messages in the buffer from being
 * deleted, but migration still dies (even on one node).
 *
 * Revision 1.1000  1997/02/06 15:59:13  ari
 * Resetting CVS to merge branches back into the main trunk.
 * We will stick to main trunk development as suggested by CVS manual.
 * We will set up tags to track fixed points of development/release
 * as suggested by CVS manual - all praise the CVS manual.
 *
 * Revision 1.778  1997/01/28 00:31:18  ari
 * internal release uplevel to 1.778
 *
 * Revision 1.777.2.1  1997/01/27 22:45:38  ari
 * Basic Atom Migration Code added.
 * Added correct magic first line to .h files for xemacs to go to C++ mode.
 * Compiles and runs without migration turned on.
 *
 * Revision 1.777  1997/01/17 19:36:53  ari
 * Internal CVS leveling release.  Start development code work
 * at 1.777.1.1.
 *
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
