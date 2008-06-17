/**
***  Copyright (c) 1995, 1996, 1997, 1998, 1999, 2000 by
***  The Board of Trustees of the University of Illinois.
***  All rights reserved.
**/

#ifndef PROXYPATCH_H
#define PROXYPATCH_H

#include "Patch.h"

class ProxyDataMsg;
class ProxyAtomsMsg;
//class ProxyAllMsg;

#define PROXYMSGNOTBUFFERED 0
#define PROXYDATAMSGBUFFERED 1
#define PROXYALLMSGBUFFERED 2

class ProxyPatch : public Patch
{
  public:

     ProxyPatch(PatchID pd);
     virtual ~ProxyPatch(void);

     void receiveAtoms(ProxyAtomsMsg*);
     void receiveData(ProxyDataMsg*);
     void receiveAll(ProxyDataMsg*);

     void setSpanningTree(int, int*, int);
     int  getSpanningTreeParent() { return parent; }
     int  getSpanningTreeChild(int *);
     inline int getSpanningTreeNChild(void) {
        #ifdef NODEAWARE_PROXY_SPANNINGTREE
            return numChild;
        #else
            return nChild; 
        #endif
     }
     ProxyCombinedResultMsg *depositCombinedResultMsg(ProxyCombinedResultMsg *);

#if CMK_PERSISTENT_COMM
  private:
     PersistentHandle localphs;
#endif
  protected:

     virtual void boxClosed(int);

  private:

     void sendResults(void);

     //"proxyMsgBufferStatus" indicates whether there's a ProxyDataMsg buffered
     // and waiting to be processed, while "curProxyMsg" points to 
     // the actual msg. This msg will be freed at the next step. --Chao Mei  
     int proxyMsgBufferStatus;
     ProxyDataMsg* curProxyMsg;
     ProxyDataMsg* prevProxyMsg;

     // for spanning tree
     ProxyCombinedResultMsg *msgCBuffer;
     int parent;
#ifdef NODEAWARE_PROXY_SPANNINGTREE
     int *children;
     int numChild;
#else
     int *child; // spanning tree for recvResults()
     int nChild;
#endif
     int nWait;
};


#endif

