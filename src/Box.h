/**
***  Copyright (c) 1995, 1996, 1997, 1998, 1999, 2000 by
***  The Board of Trustees of the University of Illinois.
***  All rights reserved.
**/

#ifndef BOX_H
#define BOX_H

#include "OwnerBox.h"

template <class Owner, class Data> class Box {

  friend class OwnerBox<Owner,Data>;

  private:

    Box(OwnerBox<Owner,Data>* o): ownerBox(o) { state = CLOSED; };

    ~Box(void) {};

    enum box_state {OPEN, CLOSED} state;
    OwnerBox<Owner,Data> *ownerBox;	

  public:

    Data* open(void) {
      if (state != OPEN) {
        state = OPEN; 
        ownerBox->openCount--;
      }
      return ownerBox->data; 
    }
    void close(Data ** const t) {
      if (state != CLOSED) {
        state = CLOSED;
        *t = NULL;

        // Trigger callback!
        if ( ! --ownerBox->closeCount ) {
	  ownerBox->close();
        }
      }
    }

};

#endif // BOX_H
