
#include <stddef.h>
#include "ckdefs.h"
#include "chare.h"
#include "c++interface.h"

template <class Owner, class Data>
OwnerBox<Owner,Data>::OwnerBox(Owner *o, void (Owner::*fn)() ) : 
    owner(o), callback(fn), numberUsers(0), 
    closeCount(0), openCount(0), data(NULL) {};


template <class Owner, class Data>
OwnerBox<Owner,Data>::~OwnerBox() {
    if (numberUsers) {
      CPrintf("OwnerBox::~OwnerBox() - still have boxes out there!\n");
    }
}
      
template <class Owner, class Data>
void OwnerBox<Owner,Data>::open(Data* d) {
      closeCount = openCount = numberUsers;
      data = d;
  }

template <class Owner, class Data>
  void OwnerBox<Owner,Data>::close() {
    if (!closeCount && !openCount) {
      data = NULL; closeCount = openCount = numberUsers;
      (owner->*callback)();
    }
    else {
      CPrintf("OwnerBox::close() - close called, but \
		closeCount %d openCount %d\n", closeCount, openCount);
    }
  }

template <class Owner, class Data>
  Box<Owner,Data>* OwnerBox<Owner,Data>::checkOut(void) {
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
  int OwnerBox<Owner,Data>::isOpen() { return (openCount); };

