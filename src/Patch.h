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
     PositionBox<Patch>* registerAvgPositionPickup(ComputeID cid, int trans = 13);
     void unregisterAvgPositionPickup(ComputeID cid,
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
     void avgPositionBoxClosed(void);

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
     PositionList  p_avg;
     Position      *positionPtr;
     Position      *avgPositionPtr;
     ForceList     f[Results::maxNumForces];
     Results	   results;
     AtomPropertiesList		a;
     AtomProperties		*atomPtr;

     PositionOwnerBox<Patch> positionBox;
     ComputeIDList              positionComputeList;
     PositionOwnerBox<Patch> avgPositionBox;
     ComputeIDList              avgPositionComputeList;
     OwnerBox<Patch,Results>    forceBox;
     ComputeIDList              forceComputeList;
     OwnerBox<Patch,AtomProperties>    atomBox;
     ComputeIDList              atomComputeList;

     virtual void boxClosed(int /* box */) = 0;
     int boxesOpen;

     void loadAtomProperties(void);
     void doGroupSizeCheck(void);

     void indexAtoms();
     int _hasNewAtoms;

  private:
  

};


#endif

