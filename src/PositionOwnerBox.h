#ifndef POSITIONOWNERBOX_H
#define POSITIONOWNERBOX_H

#include "Lattice.h"

template <class Owner>
class PositionBox;

template <class Owner>
class PositionOwnerBox {
friend class PositionBox<Owner>;
public:
  PositionOwnerBox(Owner *o, void (Owner::*fn)() );
  ~PositionOwnerBox();
      
  void open(Position* d, int n, Lattice *l );

  void close(); 

  PositionBox<Owner> *checkOut(LatticeTransform t);

  void checkIn(Box<Owner,Data> * box); 

  int isOpen() { return (openCount); };

private:
  Owner *owner;
  void (Owner::*callback)();

  Position* data;

  Position* transData[27];
  int transNeeded[27];
  Lattice *lattice;

  void callback(void);

  int numberUsers, openCount, closeCount;
};

#endif
