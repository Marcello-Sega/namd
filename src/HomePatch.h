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

#ifndef HOMEPATCH_H
#define HOMEPATCH_H

#include "ckdefs.h"
#include "chare.h"
#include "c++interface.h"

#include "NamdTypes.h"
#include "Patch.h"
#include "PatchMap.h"

#include "Sequencer.h"

#include "HomePatchTypes.h"
#include "main.h"
#include "Migration.h"

class RegisterProxyMsg;
class UnregisterProxyMsg;
class ProxyResultMsg;

class HomePatch : public Patch {
  friend PatchMgr;
  friend Sequencer;
private: // for PatchMgr to use only!!
  HomePatch(PatchID, AtomIDList, PositionList, VelocityList);
  Vector min, max;

public:

  ~HomePatch();

  void registerProxy(RegisterProxyMsg *);
  void unregisterProxy(UnregisterProxyMsg *);
  void receiveResults(ProxyResultMsg *msg);
  void useSequencer(Sequencer *sequencerPtr) {sequencer=sequencerPtr;}
  void runSequencer(int numberOfCycles = 0) { 
    sequencer->run(numberOfCycles); 
  }
  
  // methods for Sequencer to use
  void positionsReady(int doMigration=0);

  void depositMigration(PatchID, MigrationList *);
  
  void addForceToMomentum(const BigReal);
  void addVelocityToPosition(const BigReal);
  
  BigReal calcKineticEnergy();
  Vector calcMomentum();
  Vector calcAngularMomentum();
  
protected:
  virtual void boxClosed(int);
  void doAtomMigration();
  int inMigration;
  int numMlBuf;
  PatchID srcID[PatchMap::MaxOneAway];
  MigrationList *mlBuf[PatchMap::MaxOneAway];


  
private:

  // PositionList  pInit;   
  VelocityList  v; 
  // ForceList     f_short;
  // ForceList     f_long;
  
  ProxyList     proxy;
  
  Sequencer  *sequencer;

  int patchMapRead;
  void readPatchMap();

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
 *	$Author: ari $	$Locker:  $		$State: Exp $
 *	$Revision: 1.1002 $	$Date: 1997/02/11 18:51:47 $
 *
 ***************************************************************************
 * REVISION HISTORY:
 *
 * $Log: HomePatch.h,v $
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
 *
 * Revision 1.11  1996/07/15 22:57:05  gursoy
 * *** empty log message ***
 *
 * Revision 1.10  1996/07/15 21:18:01  gursoy
 * *** empty log message ***
 *
 * Revision 1.9  1996/07/10 16:19:52  gursoy
 * *** empty log message ***
 *
 * Revision 1.8  1996/07/08 21:32:38  gursoy
 * ,
 *
 * Revision 1.7  1996/06/12 16:34:46  brunner
 * *** empty log message ***
 *
 * Revision 1.6  1996/06/11 22:36:35  brunner
 * *** empty log message ***
 *
 * Revision 1.5  1996/06/11 20:07:22  gursoy
 * *** empty log message ***
 *
 * Revision 1.4  1996/06/10 22:04:14  brunner
 * *** empty log message ***
 *
 * Revision 1.3  1996/06/10 20:31:57  gursoy
 * *** empty log message ***
 *
 * Revision 1.2  1996/06/10 20:30:00  gursoy
 * *** empty log message ***
 *
 * Revision 1.1  1996/05/30 21:31:36  gursoy
 * Initial revision
 *
 ***************************************************************************/
