//-*-c++-*-
/***************************************************************************/
/*                                                                         */
/*              (C) Copyright 1996 The Board of Trustees of the            */
/*                          University of Illinois                         */
/*                           All Rights Reserved                           */
/*									   */
/***************************************************************************/

/***************************************************************************
 * DESCRIPTION:
 *
 ***************************************************************************/

#ifndef NAMDTYPES_H

#define NAMDTYPES_H

#include "Vector.h"
#include "ResizeArray.h"
#include "UniqueSortedArray.h"
#include "ResizeArrayIter.h"
#include "ResizeArrayPrimIter.h"

class Patch;
class Compute;

typedef Vector Position;
typedef Vector Velocity;
typedef Vector Force;
typedef int AtomID;
typedef int AtomType;
typedef float Mass;
typedef float Charge;

typedef double Coordinate;

struct AtomProperties
{
  AtomID id;
  AtomType type;
  Mass mass;
  Charge charge;

  // other data information
  char hydrogenGroupSize;	// 0 from group members, !0 for group parents
  char nonbondedGroupSize;	// same, but variable with strict size limits
  // Bool water;	// TRUE if water atom (O or H)  NEVER USED -JCP
  unsigned char flags;	// for fixed atoms, etc. - use with & operator

  int operator==(const AtomProperties& a) {
    return( id == a.id );
  }
};

// Definitions for AtomProperties flags
#define ATOM_FIXED	0x0001
#define GROUP_FIXED	0x0002

typedef ResizeArray<Position> PositionList;
typedef ResizeArrayIter<Position> PositionListIter;
typedef ResizeArray<Velocity> VelocityList;
typedef ResizeArrayIter<Velocity> VelocityListIter;
typedef ResizeArray<Force> ForceList;
typedef ResizeArrayIter<Force> ForceListIter;
typedef ResizeArray<AtomProperties> AtomPropertiesList;
typedef ResizeArrayIter<AtomProperties> AtomPropertiesListIter;

typedef ResizeArray<AtomID> AtomIDList;

typedef int PatchID;
typedef int ComputeID;
typedef int NodeID;

typedef ResizeArray<PatchID> PatchIDList;
typedef ResizeArray<Patch *> PatchList;

typedef UniqueSortedArray<ComputeID> ComputeIDList;
typedef ResizeArrayPrimIter<ComputeID> ComputeIDListIter;

typedef ResizeArray<Compute *> ComputeList;

// See AtomMap
struct LocalID
{
  PatchID pid;
  int index;
};

// HP compiler complains that true, false "Will be" future reserved words.
//enum Boolean
//{
//  false=0,
//  true=1
//};
#ifndef BOOLTYPE
typedef int Boolean;
#define false 0
#define true 1
#endif

#endif /* NAMDTYPES_H */

/***************************************************************************
 * RCS INFORMATION:
 *
 *	$RCSfile: NamdTypes.h,v $
 *	$Author: jim $	$Locker:  $		$State: Exp $
 *	$Revision: 1.1009 $	$Date: 1998/04/14 05:58:25 $
 *
 ***************************************************************************
 * REVISION HISTORY:
 *
 * $Log: NamdTypes.h,v $
 * Revision 1.1009  1998/04/14 05:58:25  jim
 * Added automatic correction if hgroupCutoff is too small.  No more warnings.
 * However, performance wil degrade if many groups are below cutoff size.
 *
 * Revision 1.1008  1997/12/26 23:10:53  milind
 * Made namd2 to compile, link and run under linux. Merged Templates and src
 * directoriies, and removed separate definition and declaration files for
 * templates.
 *
 * Revision 1.1007  1997/09/19 08:55:34  jim
 * Added rudimentary but relatively efficient fixed atoms.  New options
 * are fixedatoms, fixedatomsfile, and fixedatomscol (nonzero means fixed).
 * Energies will be affected, although this can be fixed with a little work.
 *
 * Revision 1.1006  1997/04/03 19:59:07  nealk
 * 1) New Fopen() which handles .Z and .gz files.
 * 2) localWaters and localNonWaters lists on each patch.
 *
 * Revision 1.1005  1997/03/19 05:50:10  jim
 * Added ComputeSphericalBC, cleaned up make dependencies.
 *
 * Revision 1.1004  1997/03/15 22:15:26  jim
 * Added ComputeCylindricalBC.  Doesn't break anything but untested and
 * cylinder is along x axis (will fix soon).
 *
 * Revision 1.1003  1997/03/04 22:37:14  ari
 * Clean up of code.  Debug statements removal, dead code removal.
 * Minor fixes, output fixes.
 * Commented some code from the top->down.  Mainly reworked Namd, Node, main.
 *
 * Revision 1.1002  1997/02/13 16:17:15  ari
 * Intermediate debuging commit - working to fix deep bug in migration?
 *
 * Revision 1.1001  1997/02/11 18:51:49  ari
 * Modified with #ifdef DPMTA to safely eliminate DPMTA codes
 * fixed non-buffering of migration msgs
 * Migration works on multiple processors
 *
 * Revision 1.1000  1997/02/06 15:58:49  ari
 * Resetting CVS to merge branches back into the main trunk.
 * We will stick to main trunk development as suggested by CVS manual.
 * We will set up tags to track fixed points of development/release
 * as suggested by CVS manual - all praise the CVS manual.
 *
 * Revision 1.779  1997/02/06 15:53:17  ari
 * Updating Revision Line, getting rid of branches
 *
 * Revision 1.778.2.1  1997/02/06 05:00:44  jim
 * Added creation of full electrostatics objects.
 *
 * Revision 1.778  1997/01/28 00:30:58  ari
 * internal release uplevel to 1.778
 *
 * Revision 1.777.2.1  1997/01/27 22:45:28  ari
 * Basic Atom Migration Code added.
 * Added correct magic first line to .h files for xemacs to go to C++ mode.
 * Compiles and runs without migration turned on.
 *
 * Revision 1.777  1997/01/17 19:36:33  ari
 * Internal CVS leveling release.  Start development code work
 * at 1.777.1.1.
 *
 * Revision 1.14  1997/01/15 17:09:43  ari
 * minor changes
 *
 * Revision 1.13  1996/12/05 01:47:07  ari
 * added == to AtomProperties definition
 *
 * Revision 1.12  1996/11/30 01:27:34  jim
 * switched to realistic ComputeType definitions
 *
 * Revision 1.11  1996/10/30 01:05:50  jim
 * added AtomPropertiesList
 *
 * Revision 1.10  1996/10/30 00:32:18  jim
 * added AtomProperties definition
 *
 * Revision 1.9  1996/10/24 18:51:09  brunner
 * Added LocalID
 *
 * Revision 1.8  1996/10/16 08:22:39  ari
 * *** empty log message ***
 *
 * Revision 1.7  1996/09/10 04:16:28  jim
 * Added iterators for all lists.
 *
 * Revision 1.6  1996/08/29 00:52:06  ari
 * *** empty log message ***
 *
 * Revision 1.5  1996/08/23 21:36:58  ari
 * *** empty log message ***
 *
 * Revision 1.4  1996/08/19 21:37:02  brunner
 * Added Coordinate
 *
 * Revision 1.3  1996/08/19 21:27:51  ari
 * .
 *
 * Revision 1.2  1996/08/16 21:42:58  brunner
 * *** empty log message ***
 *
 * Revision 1.1  1996/08/16 21:00:37  brunner
 * Initial revision
 *
 ***************************************************************************/
