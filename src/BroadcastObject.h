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

#ifndef _BCASTOBJ_H
#define _BCASTOBJ_H

#include "BroadcastMgr.h"
#include "BroadcastClient.h"

template <class T> class SimpleBroadcastObject : public BroadcastClient {

  public:

    SimpleBroadcastObject(int id=0) : BroadcastClient(id) { }
    ~SimpleBroadcastObject() { }

    T get(int tag) {
      void *buf;
      while (!(buf = (BroadcastMgr::Object())->getbuf(*this, tag))) {
        suspendFor(tag);
      }
      T tmp = *(T *)buf;
      delete (T *)buf;
      return tmp;
    }
    
    void publish(int tag,const T &t ) {
      void *buf = new T;
      *(T *)buf = t;
      BroadcastMgr::Object()->send(*this, tag, buf, sizeof(T));
    }

};

#endif
