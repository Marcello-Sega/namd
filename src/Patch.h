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

#ifndef PATCH_H
#define PATCH_H

#include "NamdTypes.h"
#include "OwnerBox.h"
#include "Box.h"
#include "PositionOwnerBox.h"
#include "PositionBox.h"
#include "UniqueSortedArray.h"
#include "Lattice.h"
#include "PatchTypes.h"

class Compute;
class Sequencer;
class PatchMap;

class LocalAtomID {
public:
    AtomID atomID;
    int index;
    LocalAtomID(AtomID a, int i) : atomID(a), index(i) {};
    LocalAtomID() {};
    ~LocalAtomID() {};
    int operator < (const LocalAtomID &a) const {
	return (atomID < a.atomID);
    }
    int operator== (const LocalAtomID &a) const {
       return (atomID == a.atomID);
    }
};

typedef UniqueSortedArray<LocalAtomID> LocalIndex ;
typedef UniqueSortedArray<int> LocalInt ;

// This the base class of homepatches and proxy patches. It maintains
// common functions of these patches. These include managing dependences
// between compute (force) objects and the patch and updating atom map.

class Patch
{
  public:

     Patch(PatchID pd);
     Patch(PatchID pd, AtomIDList al, PositionList pl);
     void loadAtoms(AtomIDList al);
     int hasNewAtoms() { return _hasNewAtoms; }
     virtual ~Patch(void) { };

     // methods for use by Compute objects
     PositionBox<Patch>* registerPositionPickup(ComputeID cid, int trans = 13);
     void unregisterPositionPickup(ComputeID cid,
				   PositionBox<Patch>**const box);
     Box<Patch,Results>* registerForceDeposit(ComputeID cid);
     void unregisterForceDeposit(ComputeID cid, Box<Patch,Results> **const box);
     Box<Patch,AtomProperties>* registerAtomPickup(ComputeID cid);
     void unregisterAtomPickup(ComputeID cid,
			       Box<Patch,AtomProperties> **const box);

     // methods for use by Sequencer or ProxyManager
     // void positionsReady(void) { positionsReady(0); }
     void positionsReady(int n=0);

     // methods for Box callbacks
     void positionBoxClosed(void);
     void forceBoxClosed(void);
     void atomBoxClosed(void);

     int getNumAtoms() { return numAtoms; }
     AtomIDList &getAtomIDList() { return (atomIDList); }

     PatchID getPatchID() { return patchID; }

     Lattice &lattice;
     Flags flags;

  protected:
     static PatchMap *patchMap;

     PatchID       patchID;
     int           numAtoms;
     AtomIDList    atomIDList;
     // LocalIndex    localIndex;  NEVER USED -JCP
     // LocalInt      localWaters;  NEVER USED -JCP
     // LocalInt      localNonWaters;  NEVER USED -JCP
     PositionList  p;
     Position      *positionPtr;
     ForceList     f[Results::maxNumForces];
     Results	   results;
     AtomPropertiesList		a;
     AtomProperties		*atomPtr;

     PositionOwnerBox<Patch> positionBox;
     ComputeIDList              positionComputeList;
     OwnerBox<Patch,Results>    forceBox;
     ComputeIDList              forceComputeList;
     OwnerBox<Patch,AtomProperties>    atomBox;
     ComputeIDList              atomComputeList;

     virtual void boxClosed(int box) { box = 0; }
     int boxesOpen;

     void loadAtomProperties(void);

     void indexAtoms();
     int _hasNewAtoms;

  private:
  

};


#endif

/***************************************************************************
 * RCS INFORMATION:
 *
 *	$RCSfile: Patch.h,v $
 *	$Author: milind $	$Locker:  $		$State: Exp $
 *	$Revision: 1.1010 $	$Date: 1997/12/26 23:10:55 $
 *
 ***************************************************************************
 * REVISION HISTORY:
 *
 * $Log: Patch.h,v $
 * Revision 1.1010  1997/12/26 23:10:55  milind
 * Made namd2 to compile, link and run under linux. Merged Templates and src
 * directoriies, and removed separate definition and declaration files for
 * templates.
 *
 * Revision 1.1009  1997/09/19 05:17:44  jim
 * Cleaned up and tweaked hydrogen-group based temporary pairlist
 * generation for roughly a 6% performance improvement.
 *
 * Revision 1.1008  1997/04/08 21:08:43  jim
 * Contant pressure now correct on multiple nodes, should work with MTS.
 *
 * Revision 1.1007  1997/04/03 19:59:12  nealk
 * 1) New Fopen() which handles .Z and .gz files.
 * 2) localWaters and localNonWaters lists on each patch.
 *
 * Revision 1.1006  1997/03/19 18:10:15  nealk
 * Added sorted hydrogen group list to molecule.
 *
 * Revision 1.1005  1997/03/12 22:06:44  jim
 * First step towards multiple force returns and multiple time stepping.
 *
 * Revision 1.1004  1997/02/28 04:47:09  jim
 * Full electrostatics now works with fulldirect on one node.
 *
 * Revision 1.1003  1997/02/26 16:53:13  ari
 * Cleaning and debuging for memory leaks.
 * Adding comments.
 * Removed some dead code due to use of Quiescense detection.
 *
 * Revision 1.1002  1997/02/11 18:51:54  ari
 * Modified with #ifdef DPMTA to safely eliminate DPMTA codes
 * fixed non-buffering of migration msgs
 * Migration works on multiple processors
 *
 * Revision 1.1001  1997/02/07 05:42:31  ari
 * Some bug fixing - atom migration on one node works
 * Atom migration on multiple nodes gets SIGSEGV
 *
 * Revision 1.1000  1997/02/06 15:59:02  ari
 * Resetting CVS to merge branches back into the main trunk.
 * We will stick to main trunk development as suggested by CVS manual.
 * We will set up tags to track fixed points of development/release
 * as suggested by CVS manual - all praise the CVS manual.
 *
 * Revision 1.779  1997/02/06 15:53:21  ari
 * Updating Revision Line, getting rid of branches
 *
 * Revision 1.778.2.1  1997/02/05 22:18:17  ari
 * Added migration code - Currently the framework is
 * there with compiling code.  This version does
 * crash shortly after migration is complete.
 * Migration appears to complete, but Patches do
 * not appear to be left in a correct state.
 *
 * Revision 1.778  1997/01/28 00:31:10  ari
 * internal release uplevel to 1.778
 *
 * Revision 1.777.2.3  1997/01/27 22:45:32  ari
 * Basic Atom Migration Code added.
 * Added correct magic first line to .h files for xemacs to go to C++ mode.
 * Compiles and runs without migration turned on.
 *
 * Revision 1.777.2.2  1997/01/24 22:00:35  jim
 * Changes for periodic boundary conditions.
 *
 * Revision 1.777.2.1  1997/01/21 23:04:47  ari
 * Basic framework for atom migration placed into code.  - Non
 * functional since it is not called.  Works currently without
 * atom migration.
 *
 * Revision 1.777  1997/01/17 19:36:46  ari
 * Internal CVS leveling release.  Start development code work
 * at 1.777.1.1.
 *
 * Revision 1.11  1996/12/12 02:58:49  jim
 * cleaned stuff up, move integration routines to HomePatch
 *
 * Revision 1.10  1996/12/11 18:45:51  nealk
 * Added stubs for addForceToMomentum and addVelocityToPosition.
 *
 * Revision 1.9  1996/12/05 01:44:16  ari
 * started toward proxy management
 *
 * Revision 1.8  1996/11/30 00:41:24  jim
 * added boxesOpen counting to support HomePatch::boxClosed()
 *
 * Revision 1.7  1996/11/01 21:20:45  ari
 * *** empty log message ***
 *
 * Revision 1.6  1996/10/30 01:16:32  jim
 * added AtomProperties structure in Patch plus boxes, passing, etc.
 *
 * Revision 1.5  1996/10/29 23:35:27  ari
 * *** empty log message ***
 *
 * Revision 1.4  1996/10/16 08:22:39  ari
 * *** empty log message ***
 *
 * Revision 1.3  1996/10/04 21:07:46  jim
 * Moved in functionality from HomePatch
 *
 * Revision 1.2  1996/09/10 03:07:14  ari
 * *** empty log message ***
 *
 * Revision 1.1  1996/08/19 22:07:49  ari
 * Initial revision
 *
 ***************************************************************************/
