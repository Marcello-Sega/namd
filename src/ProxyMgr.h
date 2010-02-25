/**
***  Copyright (c) 1995, 1996, 1997, 1998, 1999, 2000 by
***  The Board of Trustees of the University of Illinois.
***  All rights reserved.
**/

#ifndef PROXYMGR_H
#define PROXYMGR_H


#include "charm++.h"
#include "envelope.h"

#include "main.h"
#include "NamdTypes.h"
#include "PatchTypes.h"
#include "UniqueSet.h"
#include "UniqueSetIter.h"
#include "ProcessorPrivate.h"
#include "ProxyMgr.decl.h"

extern int proxySendSpanning, proxyRecvSpanning;
extern const int proxySpanDim;
extern const int inNodeProxySpanDim;

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

//1. This class represents for both msg types: one that
//is originally known as ProxyAllMsg which is sent
//at the step where atoms migrate; and the other is
//sent during the steps between two migrations.
//2. In the case of memory optimized version, the scenario
//becomes tricky as load balancer will move compute objects
//around so that new ProxyPatches will be created where
//the CompAtomExt list information is not available. If
//the step immediately after the load balancing is a normal
//step, then the CompAtomExt list info has to be resent by
//the HomePatch. Because of the current Proxy msg communication
//scheme where msg is sent to ProxyMgr first, and then retransmitted
//to ProxyPatches, there's overhead when we want to resend CompAtomExt
//list as not all the ProxyPatches that are managed by ProxyMgr are
//newly created ProxyPatches. 
//--Chao Mei
class ProxyDataMsg : public CMessage_ProxyDataMsg {
public:
  PatchID patch;
  Flags flags;

  int plLen;

  CompAtom *positionList;
  int avgPlLen;
  CompAtom *avgPositionList;
  // BEGIN LA
  int vlLen;
  CompAtom *velocityList;
  // END LA

  //1. The following field will be only
  //useful for memory optimized version.
  //2. In normal case, adding this field only
  //increases the msg length by 4 bytes which
  //can be ignored considering the current fast
  //communication network
  //--Chao Mei
  int plExtLen;
  CompAtomExt *positionExtList;

#if defined(NODEAWARE_PROXY_SPANNINGTREE) && defined(USE_NODEPATCHMGR) && (CMK_SMP) && defined(NAMDSRC_IMMQD_HACK)
  //In smp layer, the couter for msg creation and process of communication
  //thread is not included in the quiescence detection process. In addition,
  //the immediate messages from other nodes are executed on the communication
  //thread. If inside the process of immediate messages, some normal Charm++
  //messages sent out which will be processed on worker threads. Then QD will
  //be a problem that the process of the normal messages sent from communication
  //thread is recorded, but the creation of such messages (although recorded
  //in the comm thread) is virtually not recorded, i.e., not visible the 
  //QD process. So we need to artificially increase the QD counter to 
  //compensate for aforementioned msg creation loss.
  //The idea is to use the following variable to indicate the normal message
  //is sent from the communication thread inside a processing of immediate
  //message. If the variable is set, then we should increase the QD counter.
  //Chao Mei
  char isFromImmMsgCall; //hack for imm msg with QD in SMP 
#endif

  // DMK - Atom Separation (water vs. non-water)
  #if NAMD_SeparateWaters != 0
    int numWaterAtoms;  // Number of atoms in positionList (from start)
	                //   that are part of water hydrogen groups.
  #endif
#ifdef REMOVE_PROXYDATAMSG_EXTRACOPY
  //Adding padding bytes to make sure that positionList is
  //32-byte aligned which usually gives better cache performance,
  //especially on BlueGene/L machine. Otherwise, we have to
  //do the extra copy.
  //The basic method to calculate padding is to add up
  //the size of all the fields so far, including
  //the message header (the envelope) , then mod (alignment)
  // --Chao Mei
#if defined(NODEAWARE_PROXY_SPANNINGTREE) && defined(USE_NODEPATCHMGR) && (CMK_SMP) && defined(NAMDSRC_IMMQD_HACK)
 #if NAMD_SeparateWaters != 0
  char padding[(32-(sizeof(envelope)+sizeof(PatchID)+sizeof(Flags)+sizeof(isFromImmMsgCall)+4*sizeof(int)+3*sizeof(void *))%32)%32];
 #else
  char padding[(32-(sizeof(envelope)+sizeof(PatchID)+sizeof(Flags)+sizeof(isFromImmMsgCall)+3*sizeof(int)+3*sizeof(void *))%32)%32];
 #endif
#else
 #if NAMD_SeparateWaters != 0
  char padding[(32-(sizeof(envelope)+sizeof(PatchID)+sizeof(Flags)+4*sizeof(int)+3*sizeof(void *))%32)%32];
 #else
  char padding[(32-(sizeof(envelope)+sizeof(PatchID)+sizeof(Flags)+3*sizeof(int)+3*sizeof(void *))%32)%32];
 #endif
#endif

#endif

};



class ProxyResultMsg : public CMessage_ProxyResultMsg {
public:
  NodeID node;
  PatchID patch;
  ForceList forceList[Results::maxNumForces];
  static void* pack(ProxyResultMsg *msg);
  static ProxyResultMsg* unpack(void *ptr);
};

class ProxyResultVarsizeMsg: public CMessage_ProxyResultVarsizeMsg{
public:
    NodeID node;
    PatchID patch;
    int flLen[Results::maxNumForces];   

    Force *forceArr;
    //Indicate the position of the force list that has zero value
    //which is not recorded in the above force array.
    char *isZero;

    //add padding bytes to make sure the beginning 
    //of force arrays is 8-byte aligned as it is originally.
    //Therefore, we have to put the forceArr field as
    //the first variable of varsize array type
    char padding[(8-(sizeof(envelope)+sizeof(NodeID)+sizeof(PatchID)+sizeof(int)*Results::maxNumForces+2*sizeof(void *))%8)%8];   

    //The length of "fls" is Results::maxNumForces
    static ProxyResultVarsizeMsg *getANewMsg(NodeID nid, PatchID pid, int prioSize, ForceList *fls); 
};

class ProxyNodeAwareSpanningTreeMsg: public CMessage_ProxyNodeAwareSpanningTreeMsg{
public:
    PatchID patch;
    NodeID procID;
    int numNodesWithProxies;
    int *numPesOfNode;
    int *allPes;

    static ProxyNodeAwareSpanningTreeMsg *getANewMsg(PatchID pid, NodeID nid, proxyTreeNode *tree, int size);

    //For debug
    void printOut(char *tag);
};

class ProxyCombinedResultMsg : public CMessage_ProxyCombinedResultMsg {
public:
  #if defined(NODEAWARE_PROXY_SPANNINGTREE) && defined(USE_NODEPATCHMGR)
  //since this msg may be processed by comm thread in the smp mode,
  //this variable helps comm thread to find which proc will actually process it.
  NodeID destPe;
  #if CMK_SMP && defined(NAMDSRC_IMMQD_HACK)
  //Mainly for QD in the presence of the optimization of using immediate
  //message. Refer to the explanation from ProxyDataMsg for the same 
  //variable. --Chao Mei
  char isFromImmMsgCall;
  #endif
  #endif
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

class ProxyTree {       // keep track of the spanning trees
  public:
    int proxyMsgCount;
    NodeIDList *proxylist;
#ifdef NODEAWARE_PROXY_SPANNINGTREE
    //a node-aware spanning tree array, each element of which
    //is a spanning tree for all proxies of a patch
    proxyTreeNodeList *naTrees;
#else
    NodeIDList *trees;
    int *sizes;
#endif
    
  public:
    ProxyTree() {
      proxyMsgCount = 0;
      proxylist = NULL;
#ifdef NODEAWARE_PROXY_SPANNINGTREE
      naTrees = NULL;
#else
      trees = NULL;
      sizes = NULL;
#endif      
    }
    ~ProxyTree() {
    }
};

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

  void setSendSpanning();
  int  getSendSpanning();

  void setRecvSpanning();
  int  getRecvSpanning();

  void buildProxySpanningTree();
  void sendSpanningTrees();
  void sendSpanningTreeToHomePatch(int pid, int *tree, int n);
  void recvSpanningTreeOnHomePatch(int pid, int *tree, int n);
  void sendSpanningTree(ProxySpanningTreeMsg *);
  void recvSpanningTree(ProxySpanningTreeMsg *);

  void sendNodeAwareSpanningTreeToHomePatch(int pid, proxyTreeNode *tree, int n);
  void recvNodeAwareSpanningTreeOnHomePatch(ProxyNodeAwareSpanningTreeMsg *msg);
  void sendNodeAwareSpanningTree(ProxyNodeAwareSpanningTreeMsg *);
  void recvNodeAwareSpanningTree(ProxyNodeAwareSpanningTreeMsg *);
  //set the proxy patch's parent field
  void recvNodeAwareSTParent(int patch, int parent);

  void buildProxySpanningTree2();               // centralized version
  void sendProxies(int pid, int *list, int n);
  void recvProxies(int pid, int *list, int n);

#ifdef NODEAWARE_PROXY_SPANNINGTREE
  void buildNodeAwareSpanningTree0();
  static void buildSinglePatchNodeAwareSpanningTree(PatchID pid, NodeIDList &proxyList, 
                                                    proxyTreeNodeList &ptnTree, int *proxyNodeMap);
#else
  void buildSpanningTree0();
#endif

  void sendResults(ProxyResultVarsizeMsg *);
  void recvResults(ProxyResultVarsizeMsg *);
  void sendResults(ProxyResultMsg *);
  void recvResults(ProxyResultMsg *);
  void sendResults(ProxyCombinedResultMsg *);
  void recvResults(ProxyCombinedResultMsg *);
  void recvImmediateResults(ProxyCombinedResultMsg *);

  void sendProxyData(ProxyDataMsg *, int, int*);
  void recvImmediateProxyData(ProxyDataMsg *);
  void recvProxyData(ProxyDataMsg *);

  void sendProxyAll(ProxyDataMsg *, int, int*);
  void recvImmediateProxyAll(ProxyDataMsg *);
  void recvProxyAll(ProxyDataMsg *);

  static ProxyMgr *Object() { return CkpvAccess(ProxyMgr_instance); }
  
  int numProxies() { return proxySet.size(); }

  static int nodecount;
  ProxyTree &getPtree();
 
private:
  ProxySet proxySet;
  ProxyTree ptree;

  void printProxySpanningTree();
};

class NodeProxyMgr : public CBase_NodeProxyMgr
{
private:
    proxyTreeNode **proxyInfo;
    int numPatches;

    CkGroupID localProxyMgr; //a charm Group variable
    PatchMap **localPatchMaps;

public:
    NodeProxyMgr(){
        proxyInfo = NULL;
        numPatches = 0;
        localPatchMaps = new PatchMap *[CkMyNodeSize()];
    }
    ~NodeProxyMgr(){
        for(int i=0; i<numPatches; i++) {
            delete proxyInfo[i];
        }
        delete [] proxyInfo;
        delete [] localPatchMaps;
    }

    void createProxyInfo(int numPs){
        numPatches = numPs;
        proxyInfo = new proxyTreeNode *[numPs];
        memset(proxyInfo, 0, sizeof(proxyTreeNode *)*numPs);
    }
    void registerPatch(int patchID, int numPes, int *pes);
    proxyTreeNode *getPatchProxyInfo(int patchID){
        return proxyInfo[patchID];
    }

    void registerLocalProxyMgr(CkGroupID one){
        localProxyMgr = one;
    }
    const CkGroupID &getLocalProxyMgr(){
        return localProxyMgr;
    }
    void registerLocalPatchMap(int rank, PatchMap *one){
        localPatchMaps[rank] = one;
    }
    PatchMap *getLocalPatchMap(int rank){
        return localPatchMaps[rank];
    }   

    void recvImmediateProxyData(ProxyDataMsg *msg);
    void recvImmediateProxyAll(ProxyDataMsg *msg);
    void recvImmediateResults(ProxyCombinedResultMsg *);
};

#endif /* PATCHMGR_H */

