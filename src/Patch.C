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


static char ident[] = "@(#)$Header: /home/cvs/namd/cvsroot/namd2/src/Patch.C,v 1.3 1996/10/04 21:07:46 jim Exp $";

#include "ckdefs.h"
#include "chare.h"
#include "c++interface.h"

#include "Patch.h"
#include "Compute.h"

Patch::Patch(PatchID pd, AtomIDList al, PositionList pl) :
   pPtr(0), fPtr(0),
   patchID(pd), atomIDList(al), p(pl),
   pBox(this,&(Patch::pClosed)),
   fBox(this,&(Patch::fClosed))
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
   // pList.add(c);
   return pBox.checkOut();
}

void Patch::unregisterPositionPickup(ComputeID cid, Box<Patch,Position> **const box)
{
   // pList.del(c);
   pBox.checkIn(*box);
   *box = 0;
}

Box<Patch,Force>* Patch::registerForceDeposit(ComputeID cid)
{
   // fList.add(c);
   return fBox.checkOut();
}

void Patch::unregisterForceDeposit(ComputeID cid, Box<Patch,Force> **const box)
{
   // fList.del(c);
   fBox.checkIn(*box);
   *box = 0;
}

void Patch::pClosed(void)
{
   p.encap(&pPtr,numAtoms);
   this->boxClosed(0);
}

void Patch::fClosed(void)
{
   f.encap(&fPtr,numAtoms);
   this->boxClosed(1);
}

void Patch::boxClosed(int box)
{
   ;
}

void Patch::positionsReady()
{
   pPtr = p.unencap();
   pBox.open(pPtr);

   fPtr = f.unencap();
   fBox.open(fPtr);

   ComputeListIter cl(pList);
   for(cl = cl.begin(); cl != cl.end(); cl++)
   {
     (*cl)->patchReady(patchID);
   } 
}


/***************************************************************************
 * RCS INFORMATION:
 *
 *	$RCSfile: Patch.C,v $
 *	$Author: jim $	$Locker:  $		$State: Exp $
 *	$Revision: 1.3 $	$Date: 1996/10/04 21:07:46 $
 *
 ***************************************************************************
 * REVISION HISTORY:
 *
 * $Log: Patch.C,v $
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
