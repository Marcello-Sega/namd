
#include <stddef.h>
#include "ckdefs.h"
#include "chare.h"
#include "c++interface.h"

template <class Owner>
PositionOwnerBox<Owner>::OwnerBox(Owner *o, void (Owner::*fn)() ) : 
    owner(o), callback(fn), numberUsers(0), 
    closeCount(0), openCount(0), data(NULL)
{
  for ( int i = 0; i < 27; ++i )
  {
    transNeeded[i] = 0;
    transData[i] = NULL;
  }
};


template <class Owner>
PositionOwnerBox<Owner>::~OwnerBox() {
    if (numberUsers) {
      CPrintf("OwnerBox::~OwnerBox() - still have boxes out there!\n");
    }
}
      
template <class Owner>
void PositionOwnerBox<Owner>::open(Data* d, int n, Lattice *l) {
      closeCount = openCount = numberUsers;
      data = d;
      lattice = l;
      for ( int i = 0; i < 27; ++i )
	if ( transNeeded[i] ) transData[i] = lattice->create(d,n,i);
  }

template <class Owner>
  void PositionOwnerBox<Owner>::close() {
    if (!closeCount && !openCount) {
      data = NULL; closeCount = openCount = numberUsers;
      for ( int i = 0; i < 27; ++i ) lattice->destroy(&transData[i],n,i);
      (owner->*callback)();
    }
    else {
      CPrintf("OwnerBox::close() - close called, but \
		closeCount %d openCount %d\n", closeCount, openCount);
    }
  }

template <class Owner>
  Box<Owner,Position>* PositionOwnerBox<Owner>::checkOut(LatticeTransform t) {
    if (closeCount != numberUsers || openCount != numberUsers) {
      CPrintf("OwnerBox::checkOut() Tried to checkOut while in use\n");
    }
    ++numberUsers; ++closeCount; ++openCount; 
    int i = t.index();
    ++transNeeded[i];
    return (new PositionBox<Owner>(this,i));
  }

template <class Owner>
  void PositionOwnerBox<Owner>::checkIn(PositionBox<Owner> * box) {
    int i = box.index();
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
  int PositionOwnerBox<Owner>::isOpen() { return (openCount); };

