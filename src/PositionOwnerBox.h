//-*-c++-*-
#ifndef POSITIONOWNERBOX_H
#define POSITIONOWNERBOX_H

#include "NamdTypes.h"

class Lattice;

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

  PositionBox<Owner> *checkOut(int i);

  void checkIn(PositionBox<Owner> * box); 

  int isOpen() {
    return (closeCount != numberUsers || openCount != numberUsers);
  }

private:
  Owner *owner;
  void (Owner::*callback)();

  Position* data;
  int numData;

  Position* transData[27];
  int transNeeded[27];
  Lattice *lattice;

  int numberUsers, openCount, closeCount;
};

#endif


/***************************************************************************
 * RCS INFORMATION:
 *
 *	$RCSfile $
 *	$Author $	$Locker:  $		$State: Exp $
 *	$Revision: 1.1003 $	$Date: 1997/04/10 09:14:08 $
 *
 ***************************************************************************
 * REVISION HISTORY:
 *
 * $Log: PositionOwnerBox.h,v $
 * Revision 1.1003  1997/04/10 09:14:08  ari
 * Final debugging for compute migration / proxy creation for load balancing.
 * Lots of debug code added, mostly turned off now.
 * Fixed bug in PositionBox when Patch had no dependencies.
 * Eliminated use of cout and misuse of iout in numerous places.
 *                                            Ari & Jim
 *
 * Revision 1.1002  1997/03/19 11:54:51  ari
 * Add Broadcast mechanism.
 * Fixed RCS Log entries on files that did not have Log entries.
 * Added some register variables to Molecule and ComputeNonbondedExcl.C
 *
 ***************************************************************************/
