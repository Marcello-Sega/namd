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

    Box(OwnerBox<Owner,Data>* o): openCount(0), ownerBox(o), user(-1) { state = CLOSED; };
    Box(OwnerBox<Owner,Data>* o,int n): openCount(0), ownerBox(o), user(n) {
      state = CLOSED;
    };

    ~Box(void) {};

    enum box_state {OPEN, CLOSED} state;
    OwnerBox<Owner,Data> *ownerBox;	

  public:
  int user;
  int openCount;

    Data* open(void) {
      if (state != OPEN) {
        openCount++;
        state = OPEN; 
        ownerBox->openCount--;
      } else {
        //do nothing
      }
      return ownerBox->data; 
    }

    void close(Data ** const t) {
      if (state != CLOSED) {
        state = CLOSED;
        *t = NULL;

        // Trigger callback!
        --ownerBox->closeCount;
        if (ownerBox->closeCount < 0)
          CmiAbort("GBISERROR! Box closed too many times!\n");
        if ( ! ownerBox->closeCount ) {
	        ownerBox->close();
        }
      } else {
        //do nothing
      }
    }

};

#endif // BOX_H
