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

#ifndef COMPUTE_H
#define COMPUTE_H

#include "main.h"
#include "ckdefs.h"
#include "chare.h"
#include "c++interface.h"

#include "NamdTypes.h"
#include "ComputeMap.h"

class Node;
class PatchMap;

// Base class for various forms of Compute objects
// including: <linkto class=ComputeAtoms>ComputeAtoms</linkto> 
// and <linkto class=ComputePatches>ComputePatches</linkto>
class Compute {
private:
  int patchReadyCounter;
  int numPatches;

protected:
  static Node* node;
  static ComputeMap *computeMap;
  static PatchMap *patchMap;

  void enqueueWork();

public:
  const ComputeID cid;

  Compute(ComputeID c) : cid(c) { computeMap = ComputeMap::Object(); }
  virtual ~Compute() {}

  static void setNode(Node *n) { node = n; }

  void setNumPatches(int n) { patchReadyCounter = numPatches = n; }
  int getNumPatches() { return (numPatches); };

  virtual void mapReady() {};
  virtual void patchReady(void);
  virtual void patchReady(PatchID pid) { if (pid > -1) patchReady(); }
  virtual void doWork();
};

#endif
/***************************************************************************
 * RCS INFORMATION:
 *
 *	$RCSfile: Compute.h,v $
 *	$Author: ari $	$Locker:  $		$State: Exp $
 *	$Revision: 1.777 $	$Date: 1997/01/17 19:35:33 $
 *
 ***************************************************************************
 * REVISION HISTORY:
 *
 * $Log: Compute.h,v $
 * Revision 1.777  1997/01/17 19:35:33  ari
 * Internal CVS leveling release.  Start development code work
 * at 1.777.1.1.
 *
 * Revision 1.5  1996/10/29 23:36:13  ari
 * *** empty log message ***
 *
 * Revision 1.4  1996/10/22 19:16:11  ari
 * *** empty log message ***
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

