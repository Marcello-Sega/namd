
#include "PositionOwnerBox.h"

#include <stddef.h>
#include "ckdefs.h"
#include "chare.h"
#include "c++interface.h"
#include "Lattice.h"

#define MIN_DEBUG_LEVEL 4
#define DEBUGM
#include "Debug.h"

template <class Owner>
PositionOwnerBox<Owner>::PositionOwnerBox(Owner *o, void (Owner::*fn)() ) : 
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
PositionOwnerBox<Owner>::~PositionOwnerBox() {
    if (numberUsers) {
      CPrintf("PositionOwnerBox::~PositionOwnerBox() - still have boxes out there!\n");
    }
}
      
template <class Owner>
void PositionOwnerBox<Owner>::open(Position* d, int n, Lattice *l) {
      closeCount = openCount = numberUsers;
      data = d;
      lattice = l;
      for ( int i = 0; i < 27; ++i )
	if ( transNeeded[i] ) transData[i] = lattice->create(d,n,i);
      DebugM(2,transNeeded[0]<<","<<transNeeded[1]<<","<<transNeeded[2]<<","
	<<transNeeded[3]<<","<<transNeeded[4]<<","<<transNeeded[5]<<","
	<<transNeeded[6]<<","<<transNeeded[7]<<","<<transNeeded[8]<<","
	<<transNeeded[9]<<","<<transNeeded[10]<<","<<transNeeded[11]<<","
	<<transNeeded[12]<<","<<transNeeded[13]<<","<<transNeeded[14]<<","
	<<transNeeded[15]<<","<<transNeeded[16]<<","<<transNeeded[17]<<","
	<<transNeeded[18]<<","<<transNeeded[19]<<","<<transNeeded[20]<<","
	<<transNeeded[21]<<","<<transNeeded[22]<<","<<transNeeded[23]<<","
	<<transNeeded[24]<<","<<transNeeded[25]<<","<<transNeeded[26]<<"\n");
  }

template <class Owner>
  void PositionOwnerBox<Owner>::close() {
    if (!closeCount && !openCount) {
      data = NULL; closeCount = openCount = numberUsers;
      for ( int i = 0; i < 27; ++i ) lattice->destroy(&transData[i],i);
      (owner->*callback)();
    }
    else {
      CPrintf("OwnerBox::close() - close called, but \
		closeCount %d openCount %d\n", closeCount, openCount);
    }
  }

template <class Owner>
  PositionBox<Owner>* PositionOwnerBox<Owner>::checkOut(int i) {
    if (closeCount != numberUsers || openCount != numberUsers) {
      CPrintf("OwnerBox::checkOut() Tried to checkOut while in use\n");
    }
    ++numberUsers; ++closeCount; ++openCount; 
    ++transNeeded[i];
    return (new PositionBox<Owner>(this,i));
  }

template <class Owner>
  void PositionOwnerBox<Owner>::checkIn(PositionBox<Owner> * box) {
    int i = box->trans;
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


/***************************************************************************
 * RCS INFORMATION:
 *
 *	$RCSfile $
 *	$Author $	$Locker:  $		$State: Exp $
 *	$Revision: 1.1002 $	$Date: 1997/03/19 11:54:50 $
 *
 ***************************************************************************
 * REVISION HISTORY:
 *
 * $Log: PositionOwnerBox.C,v $
 * Revision 1.1002  1997/03/19 11:54:50  ari
 * Add Broadcast mechanism.
 * Fixed RCS Log entries on files that did not have Log entries.
 * Added some register variables to Molecule and ComputeNonbondedExcl.C
 *
 ***************************************************************************/
