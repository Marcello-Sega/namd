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

     // list of bonded computation objects that depends on me.
     // Currently, we have one from each bonded interaction type.
     // i.e. one bond, one angle, one dihedral and one improper computation
     // object per processing node.
     ComputeList bondedForces;
     ComputeList shortElectForces;

     // number of unfinished computation objects
     int bondedCounter;
     int shortElectCounter;

 //    void informCompObjs(ComputeList& compList);

  protected:
     int myPatchId;
     int numAtoms;

     // the global numbers of my atoms
     int atoms;  
  public :

     Patch();
     ~Patch();

     // register_compute is overloded
//     void registerCompute(BondForce *comp)   {bondedForces.add(comp)};
//     void registerCompute(AngleForce *comp)  {bondedForces.add(comp)}; 
//     void registerCompute(DihedForce *comp)  {bondedForces.add(comp)};
//     void registerCompute(ImpropForce *comp) {bondedForces.add(comp)};
//     void registerCompute(ShortElectForce *comp) {shortElectForces.add(comp);}
//     void registerCompute(FullElectForce *comp)  {fullElectForces.add(comp);}

     // compute_done is overloded 
//     void computeDone(BondForce *comp)      {bondedCounter--; fShortCheck();}

//     void computeDone(AngleForce *comp)     {bondedCounter--; fShortCheck();}

//     void computeDone(DihedForce *comp)     {bondedCounter--; fShortCheck();}

//     void computeDone(ImpropForce *comp)    {bondedCounter--; fShortCheck();}
     
//     void computeDone(FullElectForce *comp) {fullElectCounter--;fLongCheck();} 
  
//     void computeDone(ShortElectForce *comp) 
//     {
//          shortElectCounter--;
//          fShortCheck(); 
//          fLongCheck();
//     }


//     void  patchBasedRegDone();
//     void  atomBasedRegDone();

//     void  updateAtomMap();

//     Force *acquireBuffer(int size);
//     void  releaseBuffer(Force *buffer);
};


#endif

/***************************************************************************
 * RCS INFORMATION:
 *
 *	$RCSfile: Patch.h,v $
 *	$Author: ari $	$Locker:  $		$State: Exp $
 *	$Revision: 1.1 $	$Date: 1996/08/19 22:07:49 $
 *
 ***************************************************************************
 * REVISION HISTORY:
 *
 * $Log: Patch.h,v $
 * Revision 1.1  1996/08/19 22:07:49  ari
 * Initial revision
 *
 ***************************************************************************/
