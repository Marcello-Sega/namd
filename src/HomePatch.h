//-*-c++-*-
/***************************************************************************/
/*          (C) Copyright 1996, 1997 The Board of Trustees of the          */
/*                          University of Illinois                         */
/*                           All Rights Reserved                           */
/***************************************************************************/
/***************************************************************************
 * DESCRIPTION: HomePatch is the key distributed source/sink of Atom data
 *		including positions, velocities and forces applied
 *
 ***************************************************************************/
#ifndef HOMEPATCH_H
#define HOMEPATCH_H

#include "ckdefs.h"
#include "chare.h"
#include "c++interface.h"

#include "NamdTypes.h"
#include "Patch.h"
#include "PatchMap.h"

#include "MigrateAtomsMsg.h"
#include "HomePatchTypes.h"
#include "main.h"
#include "Migration.h"

class RegisterProxyMsg;
class UnregisterProxyMsg;
class ProxyResultMsg;
class Sequencer;

class HomePatch : public Patch {
  friend class PatchMgr;
  friend class Sequencer;

private: 
  static char ident[];
  // for PatchMgr to use only
  HomePatch(PatchID, AtomIDList, PositionList, VelocityList);
  ScaledPosition min, max, center;

public:
  ~HomePatch();

  // Message from ProxyPatch (via ProxyMgr) which registers its existence
  void registerProxy(RegisterProxyMsg *);
  // opposite of above
  void unregisterProxy(UnregisterProxyMsg *);

  // ProxyPatch sends Forces back to here (via ProxyMgr)
  void receiveResults(ProxyResultMsg *msg);

  // AtomMigration messages passes from neighbor HomePatches to here.
  void depositMigration(MigrateAtomsMsg *);

  // Bind a Sequencer to this HomePatch
  void useSequencer(Sequencer *sequencerPtr);
  // start simulation over this Patch of atoms
  void runSequencer(int numberOfCycles = 0);
  
  //--------------------------------------------------------------------
  // methods for Sequencer to use
  //

  // Signal HomePatch that positions stored are to be now to be used
  void positionsReady(int doMigration=0);

  // methods to implement integration
  void addForceToMomentum(const BigReal, const int ftag = Results::normal);
  void addVelocityToPosition(const BigReal);
  
  // patch-wise calculations
  BigReal calcKineticEnergy();
  Vector calcMomentum();
  Vector calcAngularMomentum();

  // load-balancing trigger
  void submitLoadStats(int timestep);

protected:
  virtual void boxClosed(int);

  // Internal Atom Migration methods and data
  void doGroupSizeCheck();
  void doMarginCheck();
  void doAtomMigration();
  int inMigration;
  int numMlBuf;
  MigrateAtomsMsg *msgbuf[PatchMap::MaxOneAway];
  
private:
  // Store of Atom-wise variables
  VelocityList  v; 
  PositionList  pInit;   
  
  // List of Proxies
  ProxyList     proxy;
  
  Sequencer  *sequencer;

  // Needed for initialization
  int patchMapRead;
  void readPatchMap();

  // Atom Migration internals
  int allMigrationIn;
  int migrationSuspended;
  int patchMigrationCounter;
  int numNeighbors;
  MigrationInfo realInfo[PatchMap::MaxOneAway];
  MigrationInfo *mInfo[3][3][3];
};

#endif
/***************************************************************************
 * RCS INFORMATION:
 *
 *	$RCSfile: HomePatch.h,v $
 *	$Author: jim $	$Locker:  $		$State: Exp $
 *	$Revision: 1.1013 $	$Date: 1998/01/15 04:58:48 $
 *
 ***************************************************************************
 * REVISION HISTORY:
 *
 * $Log: HomePatch.h,v $
 * Revision 1.1013  1998/01/15 04:58:48  jim
 * Corrected "friend foo" to "friend class foo".
 *
 * Revision 1.1012  1998/01/14 00:40:57  jim
 * Added hydrogen group size checking (vs. hgroupcutoff parameter).
 *
 * Revision 1.1011  1998/01/13 23:10:58  jim
 * Added margin checking - prelude to automatic migration.
 *
 * Revision 1.1010  1997/03/27 20:25:45  brunner
 * Changes for LdbCoordinator, the load balance control BOC
 *
 * Revision 1.1009  1997/03/27 08:04:17  jim
 * Reworked Lattice to keep center of cell fixed during rescaling.
 *
 * Revision 1.1008  1997/03/18 18:09:03  jim
 * Revamped collection system to ensure ordering and eliminate
 * unnecessary collections.  Also reduced make dependencies.
 *
 * Revision 1.1007  1997/03/12 22:06:41  jim
 * First step towards multiple force returns and multiple time stepping.
 *
 * Revision 1.1006  1997/03/10 17:40:12  ari
 * UniqueSet changes - some more commenting and cleanup
 *
 * Revision 1.1005  1997/03/06 22:06:02  ari
 * Removed Compute.ci
 * Comments added - more code cleaning
 *
 * Revision 1.1004  1997/02/26 21:39:28  jim
 * Fixed migration with periodic boundary conditions to correctly
 * re-center migrated atoms on their new home patch.
 *
 * Revision 1.1003  1997/02/17 23:47:00  ari
 * Added files for cleaning up atom migration code
 *
 * Revision 1.1002  1997/02/11 18:51:47  ari
 * Modified with #ifdef DPMTA to safely eliminate DPMTA codes
 * fixed non-buffering of migration msgs
 * Migration works on multiple processors
 *
 * Revision 1.1001  1997/02/10 08:17:30  jim
 * Commented out sending and allocation of unused data.
 *
 * Revision 1.1000  1997/02/06 15:58:27  ari
 * Resetting CVS to merge branches back into the main trunk.
 * We will stick to main trunk development as suggested by CVS manual.
 * We will set up tags to track fixed points of development/release
 * as suggested by CVS manual - all praise the CVS manual.
 *
 * Revision 1.779  1997/02/06 15:53:12  ari
 * Updating Revision Line, getting rid of branches
 *
 * Revision 1.778.2.1  1997/02/05 22:18:14  ari
 * Added migration code - Currently the framework is
 * there with compiling code.  This version does
 * crash shortly after migration is complete.
 * Migration appears to complete, but Patches do
 * not appear to be left in a correct state.
 *
 * Revision 1.778  1997/01/28 00:30:36  ari
 * internal release uplevel to 1.778
 *
 * Revision 1.777.2.2  1997/01/27 22:45:15  ari
 * Basic Atom Migration Code added.
 * Added correct magic first line to .h files for xemacs to go to C++ mode.
 * Compiles and runs without migration turned on.
 *
 * Revision 1.777.2.1  1997/01/21 23:04:45  ari
 * Basic framework for atom migration placed into code.  - Non
 * functional since it is not called.  Works currently without
 * atom migration.
 *
 * Revision 1.777  1997/01/17 19:36:10  ari
 * Internal CVS leveling release.  Start development code work
 * at 1.777.1.1.
 *
 * Revision 1.15  1997/01/15 17:09:43  ari
 * Commented out new Atom migration code
 *
 * Revision 1.14  1997/01/10 22:38:37  jim
 * kinetic energy reporting
 *
 * Revision 1.13  1996/12/17 23:58:02  jim
 * proxy result reporting is working
 *
 * Revision 1.12  1996/12/17 22:13:22  jim
 * implemented ProxyDataMsg use
 *
 * Revision 1.11  1996/12/17 17:07:41  jim
 * moved messages from main to ProxyMgr
 *
 * Revision 1.10  1996/12/11 22:31:41  jim
 * added integration methods for Sequencer
 *
 * Revision 1.9  1996/12/05 23:45:09  ari
 * *** empty log message ***
 *
 * Revision 1.8  1996/12/05 01:44:16  ari
 * started toward proxy management
 *
 * Revision 1.7  1996/12/01 22:46:11  jim
 * switched to use simParams for number of cycles
 *
 * Revision 1.6  1996/11/30 00:35:51  jim
 * implemented boxClosed(), useSequencer(), runSequencer()
 *
 * Revision 1.5  1996/10/16 08:22:39  ari
 * *** empty log message ***
 *
 * Revision 1.4  1996/10/04 21:14:47  jim
 * Moved functionality to Patch
 *
 * Revision 1.3  1996/09/03 22:54:25  ari
 * *** empty log message ***
 *
 * Revision 1.2  1996/08/29 00:50:42  ari
 * *** empty log message ***
 *
 * Revision 1.1  1996/08/19 22:07:49  ari
 * Initial revision
 ***************************************************************************/
