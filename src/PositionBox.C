
#include "PositionBox.h"

#include <stddef.h>

template <class Owner>
Position* PositionBox<Owner>::open(void) { 
    if (state != OPEN) {
      state = OPEN; 
      ownerBox->openCount--;
    }
    return ownerBox->transData[trans];
  }

// Closed access to the pointer
template <class Owner>
void PositionBox<Owner>::close(Position ** const t) {
    if (state != CLOSED) {
      state = CLOSED;
      *t = NULL;

      // Trigger callback!
      if ( ! --ownerBox->closeCount ) {
	ownerBox->close();
      }
    }
  }

template <class Owner>
PositionBox<Owner>::PositionBox(PositionOwnerBox<Owner>* o, int t)
	: ownerBox(o), trans(t) { state = CLOSED; };

template <class Owner>
PositionBox<Owner>::~PositionBox() {};
