/**
***  Copyright (c) 1995, 1996, 1997, 1998, 1999, 2000 by
***  The Board of Trustees of the University of Illinois.
***  All rights reserved.
**/

#ifndef HYDROGEN_H
#define HYDROGEN_H

#include "NamdTypes.h"
#include "SortableResizeArray.h"

// List maintaining the global atom indicies sorted by helix groups.
class HydrogenGroupID {
  public:
    AtomID atomID;      // global atom ID
    // isGP and atomsInGroup are determined when hydrogen bonds are found.
    int isGP;		// flag determining whether this atom is a group parent
    int atomsInGroup;   // positive number means parent of group.
                        // 0 means not in group.
    int waterVal;	// number of H bonded to O parent, 2 if water
    // although the Molecule object contains get_mother_atom(), we cannot
    // use it since Molecule.h, Node.h, and structure.h would have cyclical
    // include statements.
    int GPID;	// group parent ID

    HydrogenGroupID() {};
    ~HydrogenGroupID() {};

    int operator < (const HydrogenGroupID &a) const {
      int rval;
      if (isGP)	// same group.  Check for other hydrogen ordering
	{
	// case 1: both are group parents
	if (a.isGP)	rval = (atomID < a.atomID);
	// case 2: only this atom is a group parent
	else		rval = (atomID <= a.GPID);
	}
      else
	{
	// case 3: only 'a' is a group parent
	if (a.isGP)	rval = (GPID < a.atomID);
	// case 4: both are in a group
	// case 4.1: both in different groups
	else if (GPID != a.GPID)	rval = (GPID < a.GPID);
	// case 4.2: both in same group
	else rval = (atomID < a.atomID);
	}
      return(rval);
    }

    int operator == (const HydrogenGroupID &a) const {
      // only the same when both part of same group
      return (!isGP && !a.isGP && (GPID == a.GPID) );
    }
};

typedef SortableResizeArray<HydrogenGroupID> HydrogenGroup ;

#endif

