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

PatchMap *Patch::patchMap=0;

Patch::Patch(PatchID pd) :
   lattice(flags.lattice),
   patchID(pd), numAtoms(0),
   positionPtr(0), atomPtr(0),
   positionBox(this,&Patch::positionBoxClosed),
   avgPositionBox(this,&Patch::avgPositionBoxClosed),
   forceBox(this,&Patch::forceBoxClosed),
   atomBox(this,&Patch::atomBoxClosed),
   boxesOpen(0)
{
  lattice = Node::Object()->simParameters->lattice;
  _hasNewAtoms = 0;
}

Patch::Patch(PatchID pd, AtomIDList al, PositionList pl) :
   lattice(flags.lattice),
   patchID(pd), numAtoms(al.size()), atomIDList(al), p(pl),
   positionPtr(0), atomPtr(0),
   positionBox(this,&Patch::positionBoxClosed),
   avgPositionBox(this,&Patch::avgPositionBoxClosed),
   forceBox(this,&Patch::forceBoxClosed),
   atomBox(this,&Patch::atomBoxClosed),
   boxesOpen(0)
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
	 " on node " << CkMyPe() << endl);
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

PositionBox<Patch>* Patch::registerAvgPositionPickup(ComputeID cid, int trans)
{
   //DebugM(4, "registerAvgPositionPickup("<<patchID<<") from " << cid << "\n");
   return avgPositionBox.checkOut(trans);
}

void Patch::unregisterAvgPositionPickup(ComputeID cid, PositionBox<Patch> **const box)
{
   DebugM(4, "UnregisterAvgPositionPickup from " << cid << "\n");
   avgPositionBox.checkIn(*box);
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

void Patch::avgPositionBoxClosed(void)
{
   p_avg.encap(&avgPositionPtr,numAtoms);
   this->boxClosed(3);
}

// void Patch::boxClosed(int box) is virtual

void Patch::positionsReady(int doneMigration)
{
   DebugM(4,"Patch::positionsReady() - patchID(" << patchID <<")"<<endl );
   ComputeMap *computeMap = ComputeMap::Object();

   doGroupSizeCheck();

   boxesOpen = 3;
   if ( flags.doMolly ) boxesOpen++;
   _hasNewAtoms = (doneMigration != 0);

   // Give all position pickup boxes access to positions
   positionPtr = p.unencap();
   positionBox.open(positionPtr,numAtoms,&lattice);
   if ( flags.doMolly ) {
     avgPositionPtr = p_avg.unencap();
     avgPositionBox.open(avgPositionPtr,numAtoms,&lattice);
   }

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
//     iout << "patchID("<<patchID<<") has no computes dependent\n" << endi;
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

void Patch::doGroupSizeCheck(void)
{
  if ( ! flags.doNonbonded ) return;

  SimParameters *simParams = Node::Object()->simParameters;
  BigReal hgcut = 0.5 * simParams->hgroupCutoff;  hgcut *= hgcut;

  AtomPropertiesList::iterator a_i = a.begin();
  PositionList::iterator p_i = p.begin();
  PositionList::iterator p_e = p.end();

  while ( p_i != p_e ) {
    int hgs = a_i->hydrogenGroupSize;
    a_i->nonbondedGroupSize = hgs;
    ++a_i;
    BigReal x = p_i->x;
    BigReal y = p_i->y;
    BigReal z = p_i->z;
    ++p_i;
    int oversize = 0;
    for ( int i = 1; i < hgs; ++i ) {
      a_i->nonbondedGroupSize = 0;
      ++a_i;
      BigReal dx = p_i->x - x;
      BigReal dy = p_i->y - y;
      BigReal dz = p_i->z - z;
      BigReal r2 = dx * dx + dy * dy + dz * dz;
      ++p_i;
      if ( r2 > hgcut ) oversize = 1;
    }
    if ( oversize ) {
      a_i -= hgs;
      for ( int i = 0; i < hgs; ++i ) {
        a_i->nonbondedGroupSize = 1;
        ++a_i;
      }
    }
  }
}

