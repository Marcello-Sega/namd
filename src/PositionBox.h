//-*-c++-*-
#ifndef POSITIONBOX_H
#define POSITIONBOX_H

#include "NamdTypes.h"

template <class Owner>
class PositionOwnerBox;

template <class Owner>
class PositionBox {
friend class PositionOwnerBox<Owner>;
public:
  Position* open(void);
  void close(Position ** const t);

private:
  PositionBox(PositionOwnerBox<Owner>* o, int t = 13);
  ~PositionBox();

  enum box_state {OPEN, CLOSED} state;
  PositionOwnerBox<Owner> *ownerBox;
  int trans;
};

#endif // BOX_H
