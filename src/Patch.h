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
#include "Templates/OwnerBox.h"
#include "Templates/Box.h"
#include "Templates/UniqueSortedArray.h"

class Node;

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


// This the base class of homepatches and proxy patches. It maintains
// common functions of these patches. These include managing dependences
// between compute (force) objects and the patch and updating atom map.

class Patch
{
  public:

     Patch(PatchID pd, AtomIDList al, PositionList pl);
     virtual ~Patch(void) { };

     // methods for use by Compute objects
     Box<Patch,Position>* registerPositionPickup(ComputeID cid);
     void unregisterPositionPickup(ComputeID cid, Box<Patch,Position> **const box);
     Box<Patch,Force>* registerForceDeposit(ComputeID cid);
     void unregisterForceDeposit(ComputeID cid, Box<Patch,Force> **const box);
     Box<Patch,AtomProperties>* registerAtomPickup(ComputeID cid);
     void unregisterAtomPickup(ComputeID cid, Box<Patch,AtomProperties> **const box);

     // methods for use by Sequencer or ProxyManager
     void positionsReady(void);

     // methods for Box callbacks
     void positionBoxClosed(void);
     void forceBoxClosed(void);
     void atomBoxClosed(void);

     static void setNode(Node * n) { node = n; }
     int getNumAtoms() { return numAtoms; }
     const AtomIDList &getAtomIDList() { return (atomIDList); }


  protected:
     static Node* node;

     PatchID       patchID;
     int           numAtoms;
     AtomIDList    atomIDList;
     LocalIndex    localIndex;
     PositionList  p;
     Position      *positionPtr;
     ForceList     f;
     Force         *forcePtr;
     AtomPropertiesList		a;
     AtomProperties		*atomPtr;

     OwnerBox<Patch,Position> positionBox;
     ComputeIDList              positionComputeList;
     OwnerBox<Patch,Force>    forceBox;
     ComputeIDList              forceComputeList;
     OwnerBox<Patch,AtomProperties>    atomBox;
     ComputeIDList              atomComputeList;

     virtual void boxClosed(int);

  private:

};


#endif

/***************************************************************************
 * RCS INFORMATION:
 *
 *	$RCSfile: Patch.h,v $
 *	$Author: jim $	$Locker:  $		$State: Exp $
 *	$Revision: 1.6 $	$Date: 1996/10/30 01:16:32 $
 *
 ***************************************************************************
 * REVISION HISTORY:
 *
 * $Log: Patch.h,v $
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
