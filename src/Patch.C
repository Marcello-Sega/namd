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


static char ident[] = "@(#)$Header: /home/cvs/namd/cvsroot/namd2/src/Patch.C,v 1.8 1996/11/01 21:20:45 ari Exp $";

#include "ckdefs.h"
#include "chare.h"
#include "c++interface.h"

#include "Patch.h"
#include "PatchMap.h"
#include "Compute.h"

#include "ComputeMap.h"
#include "Node.h"
#include "Molecule.h"

Patch::Patch(PatchID pd, AtomIDList al, PositionList pl) :
   positionPtr(0), forcePtr(0), atomPtr(0),
   patchID(pd), atomIDList(al), p(pl),
   positionBox(this,&(Patch::positionBoxClosed)),
   forceBox(this,&(Patch::forceBoxClosed)),
   atomBox(this,&(Patch::atomBoxClosed))
{
    patchMap = PatchMap::Object();
    patchMap->registerPatch(patchID, this);

    if (atomIDList.size() != p.size())
    {
      CPrintf(
         "Patch::Patch(...) : Different numbers of Coordinates and IDs!\n");
    }
    Molecule *mol = Node::Object()->molecule;
    AtomIDListIter ai(atomIDList);
    int i = 0;
    for ( ai = ai.begin(); ai != ai.end(); ai++ )
    {
      LocalAtomID la(*ai, i++);
      localIndex.load(la);
      AtomProperties ap;
      ap.id = *ai;
      ap.type = mol->atomvdwtype(*ai);
      ap.mass = mol->atommass(*ai);
      ap.charge = mol->atomcharge(*ai);
      a.add(ap);
    }
    localIndex.sort();
    localIndex.uniq();
}


Box<Patch,Position>* Patch::registerPositionPickup(ComputeID cid)
{
   if (positionComputeList.add(cid) < 0) return NULL;
   return positionBox.checkOut();
}

void Patch::unregisterPositionPickup(ComputeID cid, Box<Patch,Position> **const box)
{
   positionComputeList.del(cid);
   positionBox.checkIn(*box);
   *box = 0;
}

Box<Patch,Force>* Patch::registerForceDeposit(ComputeID cid)
{
   if (forceComputeList.add(cid) < 0) return NULL;
   return forceBox.checkOut();
}

void Patch::unregisterForceDeposit(ComputeID cid, Box<Patch,Force> **const box)
{
   forceComputeList.del(cid);
   forceBox.checkIn(*box);
   *box = 0;
}

Box<Patch,AtomProperties>* Patch::registerAtomPickup(ComputeID cid)
{
   if (atomComputeList.add(cid) < 0) return NULL;
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
   f.encap(&forcePtr,numAtoms);
   this->boxClosed(1);
}

void Patch::atomBoxClosed(void)
{
   a.encap(&atomPtr,numAtoms);
   this->boxClosed(2);
}

void Patch::boxClosed(int box)
{
   ;
}

void Patch::positionsReady()
{
   ComputeMap *computeMap = ComputeMap::Object();

   // Give all position pickup boxes access to positions
   positionPtr = p.unencap();
   positionBox.open(positionPtr);

   // Give all force deposit boxes access to forces
   forcePtr = f.unencap();
   forceBox.open(forcePtr);

   // Give all atom properties pickup boxes access to atom properties
   atomPtr = a.unencap();
   atomBox.open(atomPtr);

   // Iterate over compute objects that need to be informed we are ready
   ComputeIDListIter cid(positionComputeList);
   for(cid = cid.begin(); cid != cid.end(); cid++)
   {
     computeMap->compute(*cid)->patchReady(patchID);
   } 
}


/***************************************************************************
 * RCS INFORMATION:
 *
 *	$RCSfile: Patch.C,v $
 *	$Author: ari $	$Locker:  $		$State: Exp $
 *	$Revision: 1.8 $	$Date: 1996/11/01 21:20:45 $
 *
 ***************************************************************************
 * REVISION HISTORY:
 *
 * $Log: Patch.C,v $
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
