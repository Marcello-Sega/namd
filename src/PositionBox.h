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


/***************************************************************************
 * RCS INFORMATION:
 *
 *	$RCSfile $
 *	$Author $	$Locker:  $		$State: Exp $
 *	$Revision: 1.1001 $	$Date: 1997/03/19 11:54:50 $
 *
 ***************************************************************************
 * REVISION HISTORY:
 *
 * $Log: PositionBox.h,v $
 * Revision 1.1001  1997/03/19 11:54:50  ari
 * Add Broadcast mechanism.
 * Fixed RCS Log entries on files that did not have Log entries.
 * Added some register variables to Molecule and ComputeNonbondedExcl.C
 *
 ***************************************************************************/
