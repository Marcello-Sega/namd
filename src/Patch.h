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

class Compute;


// This the base class of homepatches and proxy patches. It maintains
// common functions of these patches. These include managing dependences
// between compute (force) objects and the patch and updatint atom map
class Patch
{
  private:

     ComputePickupList computePickupList;
     ComputeDepositList computeDepositList;

     // number of unfinished computation objects
     int bondedCounter;
     int shortElectCounter;

     void positionsReady();

  protected:
     PatchID myPatchId;
     int numAtoms;



  public :

     Patch();
     ~Patch();

};


#endif

/***************************************************************************
 * RCS INFORMATION:
 *
 *	$RCSfile: Patch.h,v $
 *	$Author: ari $	$Locker:  $		$State: Exp $
 *	$Revision: 1.2 $	$Date: 1996/09/10 03:07:14 $
 *
 ***************************************************************************
 * REVISION HISTORY:
 *
 * $Log: Patch.h,v $
 * Revision 1.2  1996/09/10 03:07:14  ari
 * *** empty log message ***
 *
 * Revision 1.1  1996/08/19 22:07:49  ari
 * Initial revision
 *
 ***************************************************************************/
