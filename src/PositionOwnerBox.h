/**
***  Copyright (c) 1995, 1996, 1997, 1998, 1999, 2000 by
***  The Board of Trustees of the University of Illinois.
***  All rights reserved.
**/

#ifndef POSITIONOWNERBOX_H
#define POSITIONOWNERBOX_H

#include "PositionBox.h"
#include "Lattice.h"
#include "charm++.h"

template <class Owner> class PositionOwnerBox {

  friend class PositionBox<Owner>;

  public:

    PositionOwnerBox(Owner *o, void (Owner::*fn)() );

    ~PositionOwnerBox() {
      if (numberUsers) {
        CkPrintf("PositionOwnerBox::~PositionOwnerBox() - \
                   still have boxes out there!\n");
      }
    }
      
    void open(CompAtom* d, int n, Lattice *l) {
      numData = n;
      closeCount = openCount = numberUsers;
      data = d;
      lattice = l;
      for ( int i = 0; i < 27; ++i )
        if ( transNeeded[i] ) transData[i] = lattice->create(d,n,i);
      if ( ! closeCount ) close();
    }

    inline void close(void);

    inline PositionBox<Owner>* checkOut(int i);

    inline void checkIn(PositionBox<Owner> * box);

    int isOpen() {
      return (closeCount != numberUsers || openCount != numberUsers);
    }
  
  private:
    Owner *owner;
    void (Owner::*callback)();
  
    CompAtom* data;
    int numData;
  
    CompAtom* transData[27];
    int transNeeded[27];
    Lattice *lattice;
  
    int numberUsers, openCount, closeCount;
};

template <class Owner>
inline PositionOwnerBox<Owner>::PositionOwnerBox(Owner *o, void (Owner::*fn)() ) : 
  owner(o), callback(fn), data(0),
  numberUsers(0), openCount(0), closeCount(0) {
  for( int i=0; i < 27; i++ ) {
    transNeeded[i] = 0;
    transData[i] = 0;
  }
};

template <class Owner>
inline PositionBox<Owner>* PositionOwnerBox<Owner>::checkOut(int i) {
  if (closeCount != numberUsers || openCount != numberUsers) {
    CkPrintf("OwnerBox::checkOut() Tried to checkOut while in use\n");
  }
  ++numberUsers; ++closeCount; ++openCount; 
  ++transNeeded[i];
  return (new PositionBox<Owner>(this,i));
}

template <class Owner>
inline void PositionOwnerBox<Owner>::checkIn(PositionBox<Owner> * box) {
  int i = box->trans;
  delete box;
  if (closeCount != numberUsers || openCount != numberUsers) {
    CkPrintf("OwnerBox::checkIn() Tried to checkIn while in use\n");
  }
  if ( ! numberUsers-- ) {
    CkPrintf("OwnerBox::checkIn() - no registrants remaining\n");
    numberUsers = 0;
  } else {
    --transNeeded[i];
    closeCount--; openCount--;
  }
}

template <class Owner>
inline void PositionOwnerBox<Owner>::close(void) {
  if (!closeCount && !openCount) {
    data = 0; closeCount = openCount = numberUsers;
    for ( int i = 0; i < 27; ++i ) lattice->destroy(&transData[i],i);
    (owner->*callback)();
  } else {
    CkPrintf("OwnerBox::close() - close called, but \
             closeCount %d openCount %d\n", closeCount, openCount);
  }
}

#endif
