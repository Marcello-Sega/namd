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

#ifndef COMPUTEPPAIR_H
#define COMPUTEPPAIR_H

#include "main.h"
#include "ckdefs.h"
#include "chare.h"
#include "c++interface.h"

#include "NamdTypes.h"

#include "Templates/Box.h"
#include "Templates/OwnerBox.h"

class Patch;
class Node;
class PatchMap;

class ComputePatchPair : public Compute {
private:
  PatchID patchID[2];
  Patch *patch[2];
  Box<Patch,Positions> *positionBox[2];
  Box<Patch,Forces> *forceBox[2];

protected :
  int numAtoms[2];
  virtual void mapReady() { Compute::mapReady(); }
  void depositAllForces();
  virtual void doForce(Position p[2][], Force f[2][]);

public:
  ComputePatchPair(ComputeID c, PatchID pid[]);
  ~ComputePatchPair();

  virtual void doWork();
};

#endif
/***************************************************************************
 * RCS INFORMATION:
 *
 *	$RCSfile: ComputePatch.h,v $
 *	$Author: ari $	$Locker:  $		$State: Exp $
 *	$Revision: 1.1 $	$Date: 1996/10/29 22:43:35 $
 *
 ***************************************************************************
 * REVISION HISTORY:
 *
 * $Log: ComputePatch.h,v $
 * Revision 1.1  1996/10/29 22:43:35  ari
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

