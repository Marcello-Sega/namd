//-*-c++-*-
/***************************************************************************/
/*              (C) Copyright 1996,1997 The Board of Trustees of the       */
/*                          University of Illinois                         */
/*                           All Rights Reserved                           */
/***************************************************************************/
/***************************************************************************
 * DESCRIPTION: Coordinates broadcast of a data type from a Controller/Seq
 *		to all other Controller/Sequencer type objects (they must
 *		run in a thread!)
 ***************************************************************************/

#include "ckdefs.h"
#include "chare.h"
#include "c++interface.h"
#include "main.h"
#include "UniqueSet.h"
#include "UniqueSetIter.h"
#include "ProcessorPrivate.h"

#ifndef _BCASTMGR_H
#define _BCASTMGR_H

class BroadcastMsg : public comm_object {
friend class BroadcastMgr;
public:
  ~BroadcastMsg() { }
  BroadcastMsg() { msg = 0; }

  // Standard new overload for comm_object new
  void * operator new(size_t s, int i, int p) {
    return comm_object::operator new(s,i,p);
  }
  void * operator new(size_t s, int i) {return comm_object::operator new(s,i);}
  void * operator new(size_t s) { return comm_object::operator new(s); }
  void * operator new(size_t, void *ptr) { return ptr; }

  // pack and unpack functions
  void * pack (int *length) {
    *length = size + 4*sizeof(int);

    char *buffer;
    char *b = buffer = (char *)new_packbuffer(this, *length);
    memcpy(b, (void *)&size, sizeof(int)); b += sizeof(int);
    memcpy(b, (void *)&id, sizeof(int)); b += sizeof(int);
    memcpy(b, (void *)&tag, sizeof(int)); b += sizeof(int);
    memcpy(b, (void *)&node, sizeof(int)); b += sizeof(int);
    memcpy(b, msg, size); b += size;

    delete msg; // don't need [] since msg is an array of simple type
		// with no need for destructor calls
    msg = 0;

    return(buffer);
  }
    
  void unpack (void *in) {
    new((void *)this) BroadcastMsg;
    char *b = (char *)in;
    memcpy((void *)&size, b, sizeof(int)); b += sizeof(int);
    memcpy((void *)&id, b, sizeof(int)); b += sizeof(int);
    memcpy((void *)&tag, b, sizeof(int)); b += sizeof(int);
    memcpy((void *)&node, b, sizeof(int)); b += sizeof(int);
    msg = (void *)new char[size];
    memcpy((void *)msg, b, size); b += size;
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
    : tag(t), msgSize(s), counter(c), msg(m) {}
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
  BroadcastMgr(GroupInitMsg *msg) { 
    delete msg; 
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
/***************************************************************************
 * RCS INFORMATION:
 *
 *	$RCSfile $
 *	$Author $	$Locker:  $		$State: Exp $
 *	$Revision: 1.5 $	$Date: 1997/12/26 23:10:39 $
 *
 ***************************************************************************
 * REVISION HISTORY:
 *
 * $Log: BroadcastMgr.h,v $
 * Revision 1.5  1997/12/26 23:10:39  milind
 * Made namd2 to compile, link and run under linux. Merged Templates and src
 * directoriies, and removed separate definition and declaration files for
 * templates.
 *
 * Revision 1.4  1997/11/07 20:17:33  milind
 * Made NAMD to run on shared memory machines.
 *
 * Revision 1.3  1997/09/28 10:19:01  milind
 * Fixed priorities, ReductionMgr etc.
 *
 * Revision 1.2  1997/04/04 23:34:12  milind
 * Got NAMD2 to run on Origin2000.
 * Included definitions of class static variables in C files.
 * Fixed alignment bugs by using memcpy instead of assignment in
 * pack and unpack.
 *
 * Revision 1.1  1997/03/19 11:53:54  ari
 * Add Broadcast mechanism.
 * Fixed RCS Log entries on files that did not have Log entries.
 * Added some register variables to Molecule and ComputeNonbondedExcl.C
 *
 *
 ***************************************************************************/
