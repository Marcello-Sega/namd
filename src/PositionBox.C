
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

template <class Owner>
Position* PositionBox<Owner>::open(int *num) { 
    *num = ownerBox->numData;
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


/***************************************************************************
 * RCS INFORMATION:
 *
 *	$RCSfile $
 *	$Author $	$Locker:  $		$State: Exp $
 *	$Revision: 1.1002 $	$Date: 1997/04/10 09:14:06 $
 *
 ***************************************************************************
 * REVISION HISTORY:
 *
 * $Log: PositionBox.C,v $
 * Revision 1.1002  1997/04/10 09:14:06  ari
 * Final debugging for compute migration / proxy creation for load balancing.
 * Lots of debug code added, mostly turned off now.
 * Fixed bug in PositionBox when Patch had no dependencies.
 * Eliminated use of cout and misuse of iout in numerous places.
 *                                            Ari & Jim
 *
 * Revision 1.1001  1997/03/19 11:54:49  ari
 * Add Broadcast mechanism.
 * Fixed RCS Log entries on files that did not have Log entries.
 * Added some register variables to Molecule and ComputeNonbondedExcl.C
 *
 ***************************************************************************/
