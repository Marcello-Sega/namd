//-*-c++-*-
#ifndef POSITIONOWNERBOX_H
#define POSITIONOWNERBOX_H

#include "PositionBox.h"
#include "charm++.h"

class Lattice;

template <class Owner> class PositionOwnerBox {

  friend class PositionBox<Owner>;

  public:

    PositionOwnerBox(Owner *o, void (Owner::*fn)() );

    ~PositionOwnerBox() {
      if (numberUsers) {
        CPrintf("PositionOwnerBox::~PositionOwnerBox() - \
                   still have boxes out there!\n");
      }
    }
      
    void open(Position* d, int n, Lattice *l) {
      numData = n;
      closeCount = openCount = numberUsers;
      data = d;
      lattice = l;
      for ( int i = 0; i < 27; ++i )
        if ( transNeeded[i] ) transData[i] = lattice->create(d,n,i);
      if ( ! closeCount ) close();
    }

    void close(void);

    PositionBox<Owner>* checkOut(int i);

    void checkIn(PositionBox<Owner> * box);

    int isOpen() {
      return (closeCount != numberUsers || openCount != numberUsers);
    }
  
  private:
    Owner *owner;
    void (Owner::*callback)();
  
    Position* data;
    int numData;
  
    Position* transData[27];
    int transNeeded[27];
    Lattice *lattice;
  
    int numberUsers, openCount, closeCount;
};

template <class Owner>
PositionOwnerBox<Owner>::PositionOwnerBox(Owner *o, void (Owner::*fn)() ) : 
  owner(o), callback(fn), numberUsers(0), 
  closeCount(0), openCount(0), data(0) {
  for( int i=0; i < 27; i++ ) {
    transNeeded[i] = 0;
    transData[i] = 0;
  }
};

template <class Owner>
PositionBox<Owner>* PositionOwnerBox<Owner>::checkOut(int i) {
  if (closeCount != numberUsers || openCount != numberUsers) {
    CPrintf("OwnerBox::checkOut() Tried to checkOut while in use\n");
  }
  ++numberUsers; ++closeCount; ++openCount; 
  ++transNeeded[i];
  return (new PositionBox<Owner>(this,i));
}

template <class Owner>
void PositionOwnerBox<Owner>::checkIn(PositionBox<Owner> * box) {
  int i = box->trans;
  delete box;
  if (closeCount != numberUsers || openCount != numberUsers) {
    CPrintf("OwnerBox::checkIn() Tried to checkIn while in use\n");
  }
  if ( ! numberUsers-- ) {
    CPrintf("OwnerBox::checkIn() - no registrants remaining\n");
    numberUsers = 0;
  } else {
    --transNeeded[i];
    closeCount--; openCount--;
  }
}

template <class Owner>
void PositionOwnerBox<Owner>::close(void) {
  if (!closeCount && !openCount) {
    data = 0; closeCount = openCount = numberUsers;
    for ( int i = 0; i < 27; ++i ) lattice->destroy(&transData[i],i);
    (owner->*callback)();
  } else {
    CPrintf("OwnerBox::close() - close called, but \
             closeCount %d openCount %d\n", closeCount, openCount);
  }
}

#endif
