/***************************************************************************/
/*                                                                         */
/*              (C) Copyright 1996 The Board of Trustees of the            */
/*                          University of Illinois                         */
/*                           All Rights Reserved                           */
/*                                                                         */
/***************************************************************************/

/***************************************************************************
 * DESCRIPTION:
 *
 ***************************************************************************/


static char ident[] = "@(#)$Header: /home/cvs/namd/cvsroot/namd2/src/Patch.C,v 1.1019 1998/03/26 23:28:32 jim Exp $";

#include "charm++.h"

#include "Patch.h"
#include "PatchMap.h"
#include "Compute.h"

#include "ComputeMap.h"
#include "Node.h"
#include "Molecule.h"
#include "SimParameters.h"

//#define  DEBUGM
#define MIN_DEBUG_LEVEL 4
#include "Debug.h"

Patch::Patch(PatchID pd) :
   patchID(pd),
   positionPtr(0), atomPtr(0),
   positionBox(this,&Patch::positionBoxClosed),
   forceBox(this,&Patch::forceBoxClosed),
   atomBox(this,&Patch::atomBoxClosed),
   numAtoms(0), boxesOpen(0), lattice(flags.lattice)
{
  lattice = Node::Object()->simParameters->lattice;
  _hasNewAtoms = 0;
}

Patch::Patch(PatchID pd, AtomIDList al, PositionList pl) :
   patchID(pd), atomIDList(al), p(pl),
   positionPtr(0), atomPtr(0),
   positionBox(this,&Patch::positionBoxClosed),
   forceBox(this,&Patch::forceBoxClosed),
   atomBox(this,&Patch::atomBoxClosed),
   numAtoms(al.size()), boxesOpen(0), lattice(flags.lattice)
{
  lattice = Node::Object()->simParameters->lattice;
  if (atomIDList.size() != p.size()) {
    DebugM( 10, 
	    "Patch::Patch(...) : Different numbers of Coordinates and IDs!\n");
  }
  loadAtomProperties();
  _hasNewAtoms = 0;
}

void Patch::loadAtoms(AtomIDList al)
{
  atomIDList = al;
  numAtoms = atomIDList.size();
  loadAtomProperties();
  DebugM(3,"Loading " << numAtoms << " atoms into patch " << patchID <<
	 " on node " << CMyPe() << endl);
}

void Patch::loadAtomProperties(void)
{
    Molecule *mol = Node::Object()->molecule;
    int splitPatch = (Node::Object()->simParameters->splitPatch == SPLIT_PATCH_HYDROGEN);
    a.resize(0);
    // localIndex.resize(0);  NEVER USED -JCP
    register int i;
    for ( i=0; i<numAtoms; i++)
    {
      // localIndex.load(LocalAtomID(atomIDList[i], i));  NEVER USED -JCP
      AtomProperties ap;
      ap.id = atomIDList[i];
      ap.type = mol->atomvdwtype(ap.id);
      ap.mass = mol->atommass(ap.id);
      ap.charge = mol->atomcharge(ap.id);

      // make water/nonwater hydrogen group lists
      // NOTE: only group parents know the group size.
      // When no hydrogen-group migration don't use hydrogen groups.

      if ( splitPatch ) {
	if (mol->is_hydrogenGroupParent(atomIDList[i]))
	  ap.hydrogenGroupSize = mol->get_groupSize(atomIDList[i]);
	else
	  ap.hydrogenGroupSize = 0;
      }
      else ap.hydrogenGroupSize = 1;  // A group unto itself.
      // ap.water = mol->is_water(atomIDList[i]);  NEVER USED -JCP

      ap.flags = 0;
      // fixed atoms - pass one
      if ( mol->is_atom_fixed(atomIDList[i]) ) {
	ap.flags |= ATOM_FIXED;
      }

      // add the property
      a.add(ap);
    }

    int size, allfixed, j;
    for ( i=0; i<numAtoms; i+=size)
    {
      size = a[i].hydrogenGroupSize;
      allfixed = 1;
      for ( j = 0; j < size; ++j ) {
	allfixed = ( allfixed && (a[i+j].flags & ATOM_FIXED) );
      }
      if ( allfixed ) for ( j = 0; j < size; ++j ) {
	a[i+j].flags |= GROUP_FIXED;
      }
    }

    /* NEVER USED -JCP
    localIndex.sort();
    localIndex.uniq();

    // load water/nonwater lists
    localWaters.resize(0);
    localNonWaters.resize(0);
    for(i=0; i<numAtoms; i++)
    {
      if (a[i].hydrogenGroupSize) // THIS WILL NO LONGER WORK -JCP
      {
	if (a[i].water) localWaters.load(i);
	else localNonWaters.load(i);
      }
    }
    */
}

/* NEVER USED -JCP
void
Patch::indexAtoms()
{
    numAtoms = atomIDList.size();
    localIndex.resize(0);
    for ( register int i=0; i<numAtoms; i++)
    {
      localIndex.load(LocalAtomID(atomIDList[i], i));
    }
    localIndex.sort();
    localIndex.uniq();

    // load water/nonwater lists
    localWaters.resize(0);
    localNonWaters.resize(0);
    for(i=0; i<numAtoms; i++)
    {
      if (a[i].hydrogenGroupSize) // THIS WILL NO LONGER WORK -JCP
      {
	if (a[i].water) localWaters.load(i);
	else localNonWaters.load(i);
      }
    }
}
*/



PositionBox<Patch>* Patch::registerPositionPickup(ComputeID cid, int trans)
{
   //DebugM(4, "registerPositionPickupa("<<patchID<<") from " << cid << "\n");
   if (positionComputeList.add(cid) < 0)
   {
     DebugM(7, "registerPositionPickup() failed for cid " << cid << endl);
     return NULL;
   }
   return positionBox.checkOut(trans);
}

void Patch::unregisterPositionPickup(ComputeID cid, PositionBox<Patch> **const box)
{
   DebugM(4, "UnregisterPositionPickup from " << cid << "\n");
   positionComputeList.del(cid);
   positionBox.checkIn(*box);
   *box = 0;
}

Box<Patch,Results>* Patch::registerForceDeposit(ComputeID cid)
{
   if (forceComputeList.add(cid) < 0)
   {
     DebugM(7, "registerForceDeposit() failed for cid " << cid << endl);
     DebugM(7, "  size of forceCompueList " << forceComputeList.size() << endl);
     return NULL;
   }
   return forceBox.checkOut();
}

void Patch::unregisterForceDeposit(ComputeID cid, Box<Patch,Results> **const box)
{
   DebugM(4, "unregisterForceDeposit() computeID("<<cid<<")"<<endl);
   forceComputeList.del(cid);
   forceBox.checkIn(*box);
   *box = 0;
}

Box<Patch,AtomProperties>* Patch::registerAtomPickup(ComputeID cid)
{
   if (atomComputeList.add(cid) < 0)
   {
     DebugM(7, "registerAtomPickup() failed for cid " << cid << endl);
     DebugM(7, "  size of atomComputeList " << atomComputeList.size() << endl);
     return NULL;
   }
   return atomBox.checkOut();
}

void Patch::unregisterAtomPickup(ComputeID cid, Box<Patch, AtomProperties> **const box)
{
   atomComputeList.del(cid);
   atomBox.checkIn(*box);
   *box = 0;
}

void Patch::positionBoxClosed(void)
{
   p.encap(&positionPtr,numAtoms);
   this->boxClosed(0);
}

void Patch::forceBoxClosed(void)
{
   DebugM(4, "patchID("<<patchID<<") forceBoxClosed! call\n");
   for (int j = 0; j < Results::maxNumForces; ++j )
   {
     f[j].encap(&(results.f[j]),numAtoms);
   }
   this->boxClosed(1);
}

void Patch::atomBoxClosed(void)
{
   a.encap(&atomPtr,numAtoms);
   this->boxClosed(2);
}

// void Patch::boxClosed(int box) is virtual

void Patch::positionsReady(int doneMigration)
{
   DebugM(4,"Patch::positionsReady() - patchID(" << patchID <<")"<<endl );
   ComputeMap *computeMap = ComputeMap::Object();

   boxesOpen = 3;
   _hasNewAtoms = (doneMigration != 0);

   // Give all position pickup boxes access to positions
   positionPtr = p.unencap();
   positionBox.open(positionPtr,numAtoms,&lattice);

   // Give all force deposit boxes access to forces
   Force *forcePtr;
   for ( int j = 0; j < Results::maxNumForces; ++j )
   {
      f[j].resize(numAtoms);
      forcePtr = f[j].unencap();
      for(register int i=0; i<numAtoms; i++) forcePtr[i] = 0.;
      results.f[j] = forcePtr;
   }
   forceBox.open(&results);

   // Give all atom properties pickup boxes access to atom properties
   atomPtr = a.unencap();
   atomBox.open(atomPtr);

   // process computes or immediately close up boxes
   if (!positionComputeList.size()) {
   //   iout << "patchID("<<patchID<<") has no computes dependent\n" << endi;
   }
   //    positionBoxClosed();
   //    forceBoxClosed();
   //    atomBoxClosed();
   else {
     // Iterate over compute objects that need to be informed we are ready
     ComputeIDListIter cid(positionComputeList);
     for(cid = cid.begin(); cid != cid.end(); cid++)
     {
       computeMap->compute(*cid)->patchReady(patchID,doneMigration);
     } 
  }
}


/***************************************************************************
 * RCS INFORMATION:
 *
 *	$RCSfile: Patch.C,v $
 *	$Author: jim $	$Locker:  $		$State: Exp $
 *	$Revision: 1.1019 $	$Date: 1998/03/26 23:28:32 $
 *
 ***************************************************************************
 * REVISION HISTORY:
 *
 * $Log: Patch.C,v $
 * Revision 1.1019  1998/03/26 23:28:32  jim
 * Small changes for KCC port.  Altered use of strstream in ComputeFreeEnergy.
 *
 * Revision 1.1018  1998/03/03 23:05:21  brunner
 * Changed include files for new simplified Charm++ include file structure.
 *
 * Revision 1.1017  1998/01/15 18:39:40  jim
 * Tweaks for picky aCC compiler.
 *
 * Revision 1.1016  1997/09/19 08:55:35  jim
 * Added rudimentary but relatively efficient fixed atoms.  New options
 * are fixedatoms, fixedatomsfile, and fixedatomscol (nonzero means fixed).
 * Energies will be affected, although this can be fixed with a little work.
 *
 * Revision 1.1015  1997/09/19 05:17:44  jim
 * Cleaned up and tweaked hydrogen-group based temporary pairlist
 * generation for roughly a 6% performance improvement.
 *
 * Revision 1.1014  1997/05/05 15:28:26  nealk
 * Added water-water specific code to NonbondedBase.  The cutoff for the temp
 * pairlist is currently disabled.
 *
 * Revision 1.1013  1997/04/10 09:14:00  ari
 * Final debugging for compute migration / proxy creation for load balancing.
 * Lots of debug code added, mostly turned off now.
 * Fixed bug in PositionBox when Patch had no dependencies.
 * Eliminated use of cout and misuse of iout in numerous places.
 *                                            Ari & Jim
 *
 * Revision 1.1012  1997/04/08 21:08:41  jim
 * Contant pressure now correct on multiple nodes, should work with MTS.
 *
 * Revision 1.1011  1997/04/08 07:08:49  ari
 * Modification for dynamic loadbalancing - moving computes
 * Still bug in new computes or usage of proxies/homepatches.
 * Works if ldbStrategy is none as before.
 *
 * Revision 1.1010  1997/04/06 22:45:07  ari
 * Add priorities to messages.  Mods to help proxies without computes.
 * Added quick enhancement to end of list insertion of ResizeArray(s)
 *
 * Revision 1.1009  1997/04/03 19:59:11  nealk
 * 1) New Fopen() which handles .Z and .gz files.
 * 2) localWaters and localNonWaters lists on each patch.
 *
 * Revision 1.1008  1997/03/19 11:54:46  ari
 * Add Broadcast mechanism.
 * Fixed RCS Log entries on files that did not have Log entries.
 * Added some register variables to Molecule and ComputeNonbondedExcl.C
 *
 * Revision 1.1007  1997/03/12 22:06:43  jim
 * First step towards multiple force returns and multiple time stepping.
 *
 * Revision 1.1006  1997/02/11 18:51:52  ari
 * Modified with #ifdef DPMTA to safely eliminate DPMTA codes
 * fixed non-buffering of migration msgs
 * Migration works on multiple processors
 *
 * Revision 1.1005  1997/02/07 16:49:28  jim
 * Fixing bugs that affect parallel atom migration.
 *
 * Revision 1.1004  1997/02/07 05:42:30  ari
 * Some bug fixing - atom migration on one node works
 * Atom migration on multiple nodes gets SIGSEGV
 *
 * Revision 1.1003  1997/02/06 23:25:06  jim
 * Fixed bugs.
 *
 * Revision 1.1002  1997/02/06 21:20:49  jim
 * Fixed a couple of atom migration bugs.
 *
 * Revision 1.1001  1997/02/06 18:05:29  nealk
 * Modified (added some, turned off others) debug statements.
 *
 * Revision 1.1000  1997/02/06 15:59:01  ari
 * Resetting CVS to merge branches back into the main trunk.
 * We will stick to main trunk development as suggested by CVS manual.
 * We will set up tags to track fixed points of development/release
 * as suggested by CVS manual - all praise the CVS manual.
 *
 * Revision 1.779  1997/02/06 15:53:20  ari
 * Updating Revision Line, getting rid of branches
 *
 * Revision 1.778.2.2  1997/02/06 02:35:29  jim
 * Implemented periodic boundary conditions - may not work with
 * atom migration yet, but doesn't seem to alter calculation,
 * appears to work correctly when turned on.
 * NamdState chdir's to same directory as config file in argument.
 *
 * Revision 1.778.2.1  1997/02/05 22:18:16  ari
 * Added migration code - Currently the framework is
 * there with compiling code.  This version does
 * crash shortly after migration is complete.
 * Migration appears to complete, but Patches do
 * not appear to be left in a correct state.
 *
 * Revision 1.778  1997/01/28 00:31:09  ari
 * internal release uplevel to 1.778
 *
 * Revision 1.777.2.2  1997/01/24 22:00:34  jim
 * Changes for periodic boundary conditions.
 *
 * Revision 1.777.2.1  1997/01/21 23:04:47  ari
 * Basic framework for atom migration placed into code.  - Non
 * functional since it is not called.  Works currently without
 * atom migration.
 *
 * Revision 1.777  1997/01/17 19:36:45  ari
 * Internal CVS leveling release.  Start development code work
 * at 1.777.1.1.
 *
 * Revision 1.20  1997/01/15 17:09:43  ari
 * minor changes
 *
 * Revision 1.19  1996/12/19 00:37:31  jim
 * increase MIN_DEBUG_LEVEL
 *
 * Revision 1.18  1996/12/17 22:12:05  jim
 * added boxesOpen initializer
 *
 * Revision 1.17  1996/12/17 08:55:02  jim
 * added debug statement
 *
 * Revision 1.16  1996/12/12 02:58:49  jim
 * cleaned stuff up, move integration routines to HomePatch
 *
 * Revision 1.15  1996/12/11 18:45:51  nealk
 * Added stubs for addForceToMomentum and addVelocityToPosition.
 *
 * Revision 1.14  1996/12/05 01:44:16  ari
 * started toward proxy management
 *
 * Revision 1.13  1996/12/01 02:31:37  jim
 * improved debugging, fixed boxesOpen possible bug
 *
 * Revision 1.12  1996/11/30 00:41:24  jim
 * added boxesOpen counting to support HomePatch::boxClosed()
 *
 * Revision 1.11  1996/11/22 01:02:18  ari
 * *** empty log message ***
 *
 * Revision 1.10  1996/11/22 00:18:51  ari
 * *** empty log message ***
 *
 * Revision 1.9  1996/11/04 17:13:31  ari
 * *** empty log message ***
 *
 * Revision 1.8  1996/11/01 21:20:45  ari
 * *** empty log message ***
 *
 * Revision 1.7  1996/10/30 01:50:19  jim
 * atom properties list now filled on creation
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
 * Revision 1.2  1996/09/10 03:07:04  ari
 * *** empty log message ***
 *
 * Revision 1.1  1996/08/19 22:07:49  ari
 * Initial revision
 *
 * Revision 1.6  1996/07/15 22:56:48  gursoy
 * *** empty log message ***
 *
 * Revision 1.5  1996/07/15 21:17:55  gursoy
 * *** empty log message ***
 *
 * Revision 1.4  1996/07/09 21:54:34  gursoy
 * *** empty log message ***
 *
 * Revision 1.3  1996/06/24 14:14:30  gursoy
 * *** empty log message ***
 *
 * Revision 1.2  1996/06/12 16:34:23  brunner
 * *** empty log message ***
 *
 * Revision 1.1  1996/05/30 21:31:36  gursoy
 * Initial revision
 *
 ***************************************************************************/
