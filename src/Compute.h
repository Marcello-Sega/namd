//-*-c++-*-
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
  int doAtomUpdate;

protected:
  static Node* node;
  static ComputeMap *computeMap;
  static PatchMap *patchMap;

  void enqueueWork();

public:
  const ComputeID cid;

  Compute(ComputeID c) : cid(c) { 
    doAtomUpdate = false;
    computeMap = ComputeMap::Object(); 
  }
  virtual ~Compute() {}

  static void setNode(Node *n) { node = n; }

  void setNumPatches(int n) { patchReadyCounter = numPatches = n; }
  int getNumPatches() { return (numPatches); };

  virtual void initialize() {};
  virtual void atomUpdate() {};
  // virtual void patchReady(void);
  virtual void patchReady(PatchID pid) { if (pid > -1) patchReady(pid,0); }
  virtual void patchReady(PatchID, int);
  virtual int noWork(); // cleans up and returns 1 if no work to do
  virtual void doWork(); // actually does the work if noWork() returns 0
};

#endif
/***************************************************************************
 * RCS INFORMATION:
 *
 *	$RCSfile: Compute.h,v $
 *	$Author: jim $	$Locker:  $		$State: Exp $
 *	$Revision: 1.1001 $	$Date: 1997/03/12 23:59:39 $
 *
 ***************************************************************************
 * REVISION HISTORY:
 *
 * $Log: Compute.h,v $
 * Revision 1.1001  1997/03/12 23:59:39  jim
 * Added Compute::noWork() protocol to not enqueue do-nothing compute objects.
 *
 * Revision 1.1000  1997/02/06 15:57:42  ari
 * Resetting CVS to merge branches back into the main trunk.
 * We will stick to main trunk development as suggested by CVS manual.
 * We will set up tags to track fixed points of development/release
 * as suggested by CVS manual - all praise the CVS manual.
 *
 * Revision 1.779  1997/02/06 15:52:48  ari
 * Updating Revision Line, getting rid of branches
 *
 * Revision 1.778.2.1  1997/02/05 22:17:55  ari
 * Added migration code - Currently the framework is
 * there with compiling code.  This version does
 * crash shortly after migration is complete.
 * Migration appears to complete, but Patches do
 * not appear to be left in a correct state.
 *
 * Revision 1.778  1997/01/28 00:30:00  ari
 * internal release uplevel to 1.778
 *
 * Revision 1.777.2.3  1997/01/27 22:44:56  ari
 * Basic Atom Migration Code added.
 * Added correct magic first line to .h files for xemacs to go to C++ mode.
 * Compiles and runs without migration turned on.
 *
 * Revision 1.777.2.2  1997/01/27 16:52:01  nealk
 * Removed the hundreds of warnings about an unused variable.
 *
 * Revision 1.777.2.1  1997/01/21 23:04:43  ari
 * Basic framework for atom migration placed into code.  - Non
 * functional since it is not called.  Works currently without
 * atom migration.
 *
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

