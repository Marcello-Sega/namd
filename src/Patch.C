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


static char ident[] = "@(#)$Header: /home/cvs/namd/cvsroot/namd2/src/Patch.C,v 1.5 1996/10/29 23:35:27 ari Exp $";

#include "ckdefs.h"
#include "chare.h"
#include "c++interface.h"

#include "Patch.h"
#include "PatchMap.h"
#include "Compute.h"

#include "ComputeMap.h"
#include "Node.h"

Patch::Patch(PatchID pd, AtomIDList al, PositionList pl) :
   positionPtr(0), forcePtr(0),
   patchID(pd), atomIDList(al), p(pl),
   positionBox(this,&(Patch::positionBoxClosed)),
   forceBox(this,&(Patch::forceBoxClosed))
{
    if (atomIDList.size() != p.size())
    {
      CPrintf(
         "Patch::Patch(...) : Different numbers of Coordinates and IDs!\n");
    }
    AtomIDListIter a(atomIDList);
    int i = 0;
    for ( a = a.begin(); a != a.end(); a++ )
    {
      LocalAtomID la(*a, i++);
      localIndex.load(la);
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
 *	$Revision: 1.5 $	$Date: 1996/10/29 23:35:27 $
 *
 ***************************************************************************
 * REVISION HISTORY:
 *
 * $Log: Patch.C,v $
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
