/**
***  Copyright (c) 1995, 1996, 1997, 1998, 1999, 2000 by
***  The Board of Trustees of the University of Illinois.
***  All rights reserved.
**/

/*
   Coordinates broadcast of a data type from a Controller/Seq
   to all other Controller/Sequencer type objects (they must
   run in a thread!)
*/

#include "charm++.h"
#include "main.h"
#include "UniqueSet.h"
#include "UniqueSetIter.h"
#include "ProcessorPrivate.h"
#include "BroadcastMgr.decl.h"

#ifndef _BCASTMGR_H
#define _BCASTMGR_H

class BroadcastMsg : public CMessage_BroadcastMsg {
friend class BroadcastMgr;
public:
  ~BroadcastMsg() { }
  BroadcastMsg() { msg = 0; }

  // pack and unpack functions
  static void* pack(BroadcastMsg *ptr) {
    int length = ptr->size + 4*sizeof(int);
    char *buffer;
    char *b = buffer = (char *)CkAllocBuffer(ptr, length);
    memcpy(b, (void *)&(ptr->size), sizeof(int)); b += sizeof(int);
    memcpy(b, (void *)&(ptr->id), sizeof(int)); b += sizeof(int);
    memcpy(b, (void *)&(ptr->tag), sizeof(int)); b += sizeof(int);
    memcpy(b, (void *)&(ptr->node), sizeof(int)); b += sizeof(int);
    memcpy(b, ptr->msg, ptr->size); b += ptr->size;
    delete[] (char *)ptr->msg;
    delete ptr;
    return buffer;
  }
    
  static BroadcastMsg* unpack(void *ptr) {
    void *_ptr = CkAllocBuffer(ptr, sizeof(BroadcastMsg));
    BroadcastMsg *m = new (_ptr) BroadcastMsg;
    char *b = (char *)ptr;
    memcpy((void *)&(m->size), b, sizeof(int)); b += sizeof(int);
    memcpy((void *)&(m->id), b, sizeof(int)); b += sizeof(int);
    memcpy((void *)&(m->tag), b, sizeof(int)); b += sizeof(int);
    memcpy((void *)&(m->node), b, sizeof(int)); b += sizeof(int);
    m->msg = (void *)new char[m->size];
    memcpy((void *)m->msg, b, m->size); b += m->size;
    CkFreeMsg(ptr);
    return m;
  }

private:
  // Only seen by BroadcastMgr
  void *msg;
  int size;
  int id;
  int tag;
  int node;
};

class BroadcastClient;

class BroadcastClientElem {
public:
  BroadcastClientElem() {}
  BroadcastClientElem(BroadcastClient * c) : broadcastClient(c) {}
  ~BroadcastClientElem() {}

  BroadcastClient *broadcastClient;

  size_t hash() const { return (size_t)broadcastClient; }
  int operator==(const BroadcastClientElem &b) const { 
    return broadcastClient == b.broadcastClient; 
  }
};

class TaggedMsg {
public:
  TaggedMsg() {}
  TaggedMsg(int t) : tag(t) {}
  TaggedMsg(int t, int s, int c, void *m) 
    : tag(t), counter(c), msg(m), msgSize(s) {}
  ~TaggedMsg() {}

  int tag;
  int counter;
  void *msg;
  int msgSize;

  int hash() const { return tag; }
  int operator==(const TaggedMsg &tm) const { return(tag == tm.tag); }
};

class BOID {
public:
  BOID() {}
  BOID(int id) { this->id = id; }
  ~BOID() {}

  int hash() const { return id; }
  int operator==(const BOID &b) const { return id == b.id; }
  int id;

  UniqueSet<BroadcastClientElem> *broadcastSet;
  UniqueSet<TaggedMsg> *taggedMsg;
};

class BroadcastMgr : public BOCclass
{
public:
  BroadcastMgr() { 
    CpvAccess(BroadcastMgr_instance) = this; 
  }
  ~BroadcastMgr(void);
	  
  // Singleton Access method
  inline static BroadcastMgr *Object() {
    return CpvAccess(BroadcastMgr_instance);
  }

  void *getbuf(BroadcastClient &b, int tag);
  void send(BroadcastClient &b, int tag, void *buf, size_t);
  void subscribe(BroadcastClient &bc);
  void unsubscribe(BroadcastClient &bc);
  void recvBroadcast(BroadcastMsg *msg);

private:
  UniqueSet<BOID> boid;
};

#endif /* _BCASTMGR_H */

