#ifndef OWNERBOX_H
#define OWNERBOX_H

#include "ckdefs.h"

template <class Owner, class Data> class Box;

template <class Owner, class Data> class OwnerBox {

  friend class Box<Owner,Data>;

  public:

    OwnerBox(Owner *o, void (Owner::*fn)() ) :
      owner(o), callback(fn), numberUsers(0), 
      closeCount(0), openCount(0), data(0) {};

    ~OwnerBox(void) {
      if (numberUsers) {
        CPrintf("OwnerBox::~OwnerBox() - still have boxes out there!\n");
      }
    }
        
    void open(Data* d) {
      closeCount = openCount = numberUsers;
      data = d;
      if ( ! closeCount ) close();
    }
  
    void close(void);

    Box<Owner,Data> *checkOut(void);

    void checkIn(Box<Owner,Data> * box);
  
    int isOpen() {
      return (closeCount != numberUsers || openCount != numberUsers);
    }

  private:
    Owner *owner;
    void (Owner::*callback)(void);
    Data* data;
    int numberUsers, openCount, closeCount;
};

template <class Owner, class Data>
Box<Owner,Data> *OwnerBox<Owner,Data>::checkOut(void) {
  if (closeCount != numberUsers || openCount != numberUsers) {
    CPrintf("OwnerBox::checkOut() Tried to checkOut while in use\n");
  }
  ++numberUsers; ++closeCount; ++openCount; 
  return (new Box<Owner,Data>(this));
}

template <class Owner, class Data>
void OwnerBox<Owner,Data>::checkIn(Box<Owner,Data> * box) {
  delete box;
  if (closeCount != numberUsers || openCount != numberUsers) {
    CPrintf("OwnerBox::checkIn() Tried to checkIn while in use\n");
  }
  if ( ! numberUsers-- ) {
    CPrintf("OwnerBox::checkIn() - no registrants remaining\n");
    numberUsers = 0;
  } else {
    closeCount--; openCount--;
  }
}

template <class Owner, class Data>
void OwnerBox<Owner,Data>::close(void) {
  if (!closeCount && !openCount) {
    data = 0; closeCount = openCount = numberUsers;
    (owner->*callback)();
  } else {
    CPrintf("OwnerBox::close() - close called, but \
             closeCount %d openCount %d\n", closeCount, openCount);
  }
}

#endif
