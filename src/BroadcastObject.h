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

#ifndef _BCASTOBJ_H
#define _BCASTOBJ_H

#include "BroadcastMgr.h"
#include "BroadcastClient.h"
#include "common.h"

template <class T> class SimpleBroadcastObject : public BroadcastClient {

  public:

    SimpleBroadcastObject(int id) : BroadcastClient(id) { 
      if ( sizeof(T) > BCASTMSGSIZE ) {
        NAMD_bug("SimpleBroadcastObject instantiated on class larger than BCASTMSGSIZE");
      }
    }
    ~SimpleBroadcastObject() { }

    T get(int tag) {
      T tmp;
      while ( BroadcastMgr::Object()->getbuf(*this, tag, (void*)(&tmp)) < 0 ) {
        suspendFor(tag);
      }
      return tmp;
    }
    
    void publish(int tag,const T &t ) {
      BroadcastMgr::Object()->send(*this, tag, (void*)(&t), sizeof(T));
    }

};

#endif
