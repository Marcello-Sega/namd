/**
***  Copyright (c) 1995, 1996, 1997, 1998, 1999, 2000 by
***  The Board of Trustees of the University of Illinois.
***  All rights reserved.
**/

/*
   ObjectArena template - semi-contiguous storage management
*/

#ifndef OBJECTARENA_H
#define OBJECTARENA_H

#include "ResizeArray.h"

template <class Type> class ObjectArena {
  public:
    ObjectArena(void) : blockSize(1024), pos(0), end(0) {};
    ~ObjectArena(void) {
      int i;
      for( i = 0; i < blocks.size(); ++i ) delete [] blocks[i];
    }
    void setBlockSize(int n) { blockSize = n; }
    inline Type* getNewArray(int n) {
      Type *rpos = pos;
      if ( n > (blockSize/2) ) {
        rpos = new Type[n];
        blocks.add(rpos);
      } else if ( ( pos += n ) > end ) {
        rpos = pos = new Type[blockSize];
        blocks.add(pos);
        end = rpos + blockSize;
        pos += n;
      }
      return rpos;
    }

  private:
    int blockSize;
    ResizeArray<Type*> blocks;
    Type *pos, *end;
};

#endif
