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

#ifndef COMPUTEANGLE_H
#define COMPUTEANGLE_H

/*

from class Parameters
void get_angle_params(Real *k, Real *theta0, Real *k_ub, Real *r_ub,
			      Index index)
{
  *k = angle_array[index].k;
  *theta0 = angle_array[index].theta0;
  *k_ub = angle_array[index].k_ub;
  *r_ub = angle_array[index].r_ub;
}

from structures.h

typedef struct angle
{
    int atom1;
    int atom2;
    int atom3;
    Index angle_type;
} Angle;

from Molecule.h
    LintList *get_angles_for_atom(int anum)
		{return (&(anglesByAtom[anum]));}


   #include "LintList.h"
   method in LintList
   int head()
   check for LIST_EMPTY
   int next()
   returns integer value stored
*/


#include "main.h"
#include "ckdefs.h"
#include "chare.h"
#include "c++interface.h"

#include "NamdTypes.h"
#include "common.h"
#include "structures.h"
#include "Compute.h"
#include "Patch.h"

#include "Templates/Box.h"
#include "Templates/OwnerBox.h"

class AnglePatchElem;

class AngleElem {
public:
    AtomID atomID[3];
    int    localIndex[3];
    AnglePatchElem *p[3];
    Index angleType;


  AngleElem() {
    atomID[0] = -1;
    atomID[1] = -1;
    atomID[2] = -1;
    p[0] = NULL;
    p[1] = NULL;
    p[2] = NULL;
  }
  AngleElem(const Angle *a) {
    atomID[0] = a->atom1;
    atomID[1] = a->atom2;
    atomID[2] = a->atom3;
    angleType = a->angle_type;
  }

  AngleElem(AtomID atom0, AtomID atom1, AtomID atom2) {
    if (atom0 > atom2) {  // Swap end atoms so lowest is first!
      AtomID tmp = atom2; atom2 = atom0; atom0 = tmp; 
    }
    atomID[0] = atom0;
    atomID[1] = atom1;
    atomID[2] = atom2;
  }
  ~AngleElem() {};

  int operator==(const AngleElem &a) const {
    return (a.atomID[0] == atomID[0] && a.atomID[1] == atomID[1] &&
        a.atomID[2] == atomID[2]);
  }

  int operator<(const AngleElem &a) const {
    return (atomID[0] < a.atomID[0] ||
            (atomID[0] == a.atomID[0] &&
            (atomID[1] < a.atomID[1] ||
            (atomID[1] == a.atomID[1] &&
             atomID[2] < a.atomID[2]) )));
  }
};

typedef UniqueSortedArray<AngleElem> AngleList;

enum PatchType {HOME,PROXY};

class AnglePatchElem {
  public:
    PatchID patchID;
    Patch *p;
    PatchType patchType;
    Box<Patch,Position> *positionBox;
    Box<Patch,Force> *forceBox;
    Position *x;
    Force *f;

  AnglePatchElem() {
    patchID = -1;
    p = NULL;
    positionBox = NULL;
    forceBox = NULL;
    x = NULL;
    f = NULL;
  }

  AnglePatchElem(PatchID p) {
    patchID = p;
  }

  AnglePatchElem(Patch *p, PatchType pt, ComputeID cid) {
    patchID = p->getPatchID();
    this->p = p;
    patchType = pt;
    positionBox = p->registerPositionPickup(cid);
    forceBox = p->registerForceDeposit(cid);
    x = NULL;
    f = NULL;
  }
    
  ~AnglePatchElem() {};

  int operator==(const AnglePatchElem &a) const {
    return (a.patchID == patchID);
  }

  int operator<(const AnglePatchElem &a) const {
    return (patchID < a.patchID);
  }
};

typedef UniqueSortedArray<AnglePatchElem> AnglePatchList;

class AtomMap;

class ComputeAngles : public Compute {
private:
  AngleList angleList;
  AnglePatchList anglePatchList;

  PatchMap *patchMap;
  AtomMap *atomMap;

  int maxProxyAtoms;
  Force *dummy;
  
  BigReal angleForce(const Position p1, const Position p2, const Position p3,
		  Force *f1, Force *f2, Force *f3,
		    const Index angleType);

public:
  ComputeAngles(ComputeID c);
  virtual ~ComputeAngles() {
    delete [] dummy;
  }

  void mapAtoms();
  void doWork();
};

#endif
/***************************************************************************
 * RCS INFORMATION:
 *
 *	$RCSfile: ComputeHomeTuples.h,v $
 *	$Author: nealk $	$Locker:  $		$State: Exp $
 *	$Revision: 1.4 $	$Date: 1996/11/04 20:06:17 $
 *
 ***************************************************************************
 * REVISION HISTORY:
 *
 * $Log: ComputeHomeTuples.h,v $
 * Revision 1.4  1996/11/04 20:06:17  nealk
 * Now it compiles :-)
 *
 * Revision 1.3  1996/11/04 19:29:02  nealk
 * Added angleForce() to system, but it is untested.
 *
 * Revision 1.2  1996/11/04 16:55:46  ari
 * *** empty log message ***
 *
 * Revision 1.1  1996/11/01 21:20:45  ari
 * Initial revision
 *
 * Revision 1.3  1996/10/16 08:22:39  ari
 * *** empty log message ***
 *
 * Revision 1.1  1996/08/19 22:07:49  ari
 * Initial revision
 *
 * Revision 1.4  1996/07/16 01:54:12  ari
 * *** empty log message ***
 *
 * Revision 1.3  96/07/16  01:10:26  01:10:26  ari (Aritomo Shinozaki)
 * Fixed comments, added methods
 * 
 * Revision 1.2  1996/06/25 21:10:48  gursoy
 * *** empty log message ***
 *
 * Revision 1.1  1996/06/24 14:12:26  gursoy
 * Initial revision
 *
 ***************************************************************************/

