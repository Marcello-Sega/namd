//-*-c++-*-
#ifndef POSITIONBOX_H
#define POSITIONBOX_H

#include "NamdTypes.h"

template <class Owner> class PositionOwnerBox;

template <class Owner> class PositionBox {

  friend class PositionOwnerBox<Owner>;

  private:

    PositionBox(PositionOwnerBox<Owner>* o, int t=13) : ownerBox(o), trans(t) 
      { state = CLOSED; }

    ~PositionBox() {}
  
    enum box_state {OPEN, CLOSED} state;
    PositionOwnerBox<Owner> *ownerBox;
    int trans;

  public:

    Position* open(void) { 
      if (state != OPEN) {
        state = OPEN; 
        ownerBox->openCount--;
      }
      return ownerBox->transData[trans];
    }

    Position* open(int *num) { 
      *num = ownerBox->numData;
      if (state != OPEN) {
        state = OPEN; 
        ownerBox->openCount--;
      }
      return ownerBox->transData[trans];
    }

    // Closed access to the pointer
    void close(Position ** const t) {
      if (state != CLOSED) {
        state = CLOSED;
        *t = 0;
    
        // Trigger callback!
        if ( ! --ownerBox->closeCount ) {
          ownerBox->close();
         }
      }
    }

};

#endif // BOX_H
