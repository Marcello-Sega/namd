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


static char ident[] = "@(#)$Header: /home/cvs/namd/cvsroot/namd2/src/Patch.C,v 1.1 1996/08/19 22:07:49 ari Exp $";


#include "Patch.h"

//#include "Compute.h"
//#include "BondForce.h"
//#include "AngleForce.h"
//#include "DihedForce.h"
//#include "ImpropForce.h"
//#include "ElectForce.h"

Patch::Patch()
{
}

Patch::~Patch()
{
}

/*
void Patch::fShortCheck()
{
    if ( atomBasedDone      && 
         patchBasedDone     &&
         bondedCounter==0   && 
         shortElectCounter==0 ) 
    {
       bondedCounter     = numBonded;
       shortElectCounter = numShortElect;
       fShortDone();
    }
}
*/


/*
void Patch::fLongCheck()
{
    if ( patchBasedDone && shortElectCounter==0 && fullElectCounter==0)
    {
       shortElectCounter  = numShortElect;
       fullElectCounterm  = numFullElect;
       fLongDone();
    }
}
*/



/*
void Patch::patchBasedRegDone()
{
     patchBasedDone = TRUE;

     numShorElect = shortElectForces->size();
     numFullElect = fullElectForces->size();

     if ( bondedCounter==0     && shortElectCounter==0 ) fShortDone();
     if ( shortElectCounter==0 && fullElectCounter==0  ) fLongDone();
}
*/


/*
void Patch::atomBasedRegDone()
{
     atomBasedDone = TRUE;

     numBonded = bondedForces->size();

     if ( bondedCounter==0     && shortElectCounter==0 ) fShortDone();
}
*/


   /*
//void Patch::informCompObjs(ComputeList& compList)
//{
   int i;
   for(i=0; i<compList->size(); i++)
   {
      (compList[i])->patchReady();
   } 
//}
   */


/*
//void Patch::updateAtomMap()
//{
   int i;
   for(i=0; i<numAtoms;i++) 
   {
       atomMap->setAtomEntry(atoms[i],myPatchId,i);
   }
//}
*/



//Force *Patch::acquireBuffer(int size);
//void  Patch::releaseBuffer(Force *buffer);

/***************************************************************************
 * RCS INFORMATION:
 *
 *	$RCSfile: Patch.C,v $
 *	$Author: ari $	$Locker:  $		$State: Exp $
 *	$Revision: 1.1 $	$Date: 1996/08/19 22:07:49 $
 *
 ***************************************************************************
 * REVISION HISTORY:
 *
 * $Log: Patch.C,v $
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
