
#ifndef HYDROGEN_H
#define HYDROGEN_H

#include "NamdTypes.h"
#include "Templates/UniqueSortedArray.h"

// List maintaining the global atom indicies sorted by helix groups.
class HydrogenGroupID {
  public:
    AtomID atomID;      // global atom ID
    // isGP and atomsInGroup are determined when hydrogen bonds are found.
    int isGP;		// flag determining whether this atom is a group parent
    int atomsInGroup;   // positive number means parent of group.
                        // 0 means not in group.
    // although the Molecule object contains get_mother_atom(), we cannot
    // use it since Molecule.h, Node.h, and structure.h would have cyclical
    // include statements.
    int GPID;	// group parent ID

    HydrogenGroupID() {};
    ~HydrogenGroupID() {};

    int operator < (const HydrogenGroupID &a) const {
      if (isGP)
	{
	// case 1: both are group parents
	if (a.isGP)	return(atomID < a.atomID);
	// case 2: only this atom is a group parent
        return(GPID < a.atomID);
	}
      else
	{
	// case 3: only 'a' is a group parent
	if (a.isGP)	return(GPID < a.atomID);
	// case 4: both are in a group
        return(atomID < a.atomID);
	}
    }

    int operator == (const HydrogenGroupID &a) const {
      // only the same when both part of same group
      if (!isGP && !a.isGP
	  && (GPID == a.GPID) )
		return(1);	// they are the same
      return(0);
    }
};

typedef UniqueSortedArray<HydrogenGroupID> HydrogenGroup ;

#endif

/***************************************************************************
 * RCS INFORMATION:
 *
 *      $RCSfile: Hydrogen.h,v $
 *      $Author: nealk $        $Locker:  $             $State: Exp $
 *      $Revision: 1.2 $     $Date: 1997/03/19 18:47:30 $
 *
 ***************************************************************************
 * REVISION HISTORY:
 *
 * $Log: Hydrogen.h,v $
 * Revision 1.2  1997/03/19 18:47:30  nealk
 * Added log info to Hydrogen.h
 * Fixed ComputeDPMTA.C so node 0 initializes before any other nodes register
 * with the DPMTA library.
 *
 * Revision 1.1001  1997/03/19 18:10:15  nealk
 * Added sorted hydrogen group list to molecule.
 *
 ***************************************************************************

