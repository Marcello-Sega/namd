#ifndef BOX_H
#define BOX_H

#include "ckdefs.h"
#include "char.h"
#include "c++interface.h"

template <class Owner, class Data>
class OwnerBox {
friend Box<Owner,Data>
public:
  OwnerBox(Owner *o, void (Owner::*fn)()) : owner(o), callback(fn),
     numberUsers(0), numberUnclosed(0), data(NULL) {};
  ~OwnerBox() {
    if (numberUsers) {
      CPrintf("OwnerBox::~OwnerBox() - still have boxes out there!\n");
    }
  }
      
  void open(Data* d, void *box_id = 0) {
      closeCount = openCount = numberUsers;
      data = d;
  }

  void close() {
    if (!closeCount && !openCount) {
      data = NULL; closeCount = openCount = numberUsers;

    }
    else {
      CPrintf("OwnerBox::close() - close called, but \
		closeCount %d openCount %d\n", closeCount, openCount);
    }
  }

  Box<Owner,Data> *register(void) {
    if (closeCount != numberUsers || openCount != numberUsers) {
      CPrintf("OwnerBox::register() Tried to register while in use\n");
    }
    ++numberUsers; ++closeCount; ++openCount; 
    return (new Box<Owner,Data>(this));
  }

  void unregister(Box<Owner,Data> * box) {
    delete box;
    if (closeCount != numberUsers || openCount != numberUsers) {
      CPrintf("OwnerBox::unregister() Tried to unregister while in use\n");
    }
    if ( ! numberUsers-- ) {
      CPrintf("OwnerBox::unregister() - no registrants remaining\n");
      numberUsers = 0;
    } else {
      closeCount--; openCount--;
    }
  }

  int isOpen() { return (openCount); };

private:
  Owner *owner;
  void (Owner::*callback)();
  int numberUsers;
  int openCount, closeCount;
  T* data;
};


template <class Owner, class Data>
class Box {
friend OwnerBox<Owner,Data>
public:

  // Get access to a pointer
  T* open(void) { 
    if (state != OPEN) {
      state = OPEN; 
      ownerBox->openCount--;
    }
    return ownerBox->data; 
  }

  // Closed access to the pointer
  void close(void) {
    if (state != CLOSED) {
      state = CLOSED;

      // Trigger callback!
      if ( ! --ownerBox->closeCount ) {
	ownerBox->close();
      }
    }
  }

private:
  Box(OwnerBox<Owner,Data>* o) : ownerBox(o) { state = CLOSED; };
  ~Box() {};

  enum box_state {OPEN, CLOSED} state;
  OwnerBox<Owner,Data> *ownerBox;	
};

#endif // BOX_H
