/***************************************************************************/
/*                                                                         */
/*              (C) Copyright 1996 The Board of Trustees of the            */
/*                          University of Illinois                         */
/*                           All Rights Reserved                           */
/*                                                                         */
/***************************************************************************/

/***************************************************************************
 * DESCRIPTION: HomePatch owns the actual atoms of a Patch of space
 *		Proxy(s) get messages via ProxyMgr from HomePatch(es)
 *		to update lists of atoms and their coordinates
 *              HomePatch(es) also have a Sequencer bound to them
 *
 * superclass: 	Patch		
 ***************************************************************************/


#include "ckdefs.h"
#include "chare.h"
#include "c++interface.h"
#include "SimParameters.h"
#include "HomePatch.h"
#include "AtomMap.h"
#include "Node.h"
#include "PatchMap.inl"
#include "main.h"
#include "ProxyMgr.top.h"
#include "ProxyMgr.h"
#include "Migration.h"
#include "Molecule.h"
#include "PatchMgr.h"
#include "Sequencer.h"
#include "LdbCoordinator.h"

#define MIN_DEBUG_LEVEL 4
//#define DEBUGM
#include "Debug.h"

// avoid dissappearence of ident?
char HomePatch::ident[] = "@(#)$Header: /home/cvs/namd/cvsroot/namd2/src/HomePatch.C,v 1.1037 1997/10/06 00:12:31 jim Exp $";

HomePatch::HomePatch(PatchID pd, AtomIDList al, PositionList pl, 
		     VelocityList vl) : Patch(pd,al,pl), v(vl) 
{ 
  DebugM(4, "HomePatch("<<pd<<") at " << this << "\n");
  if (atomIDList.size() != v.size()) {
    CPrintf("HomePatch::HomePatch(...) : size mismatch-Velocities and IDs!\n");
  }
  AtomMap::Object()->registerIDs(pd,al);  
  min.x = PatchMap::Object()->minX(patchID);
  min.y = PatchMap::Object()->minY(patchID);
  min.z = PatchMap::Object()->minZ(patchID);
  max.x = PatchMap::Object()->maxX(patchID);
  max.y = PatchMap::Object()->maxY(patchID);
  max.z = PatchMap::Object()->maxZ(patchID);
  center = 0.5*(min+max);

  migrationSuspended = false;
  allMigrationIn = false;
  patchMapRead = 0; // We delay read of PatchMap data
		    // to make sure it is really valid
  inMigration = false;
  numMlBuf = 0;
}

// Bind a Sequencer to this HomePatch
void HomePatch::useSequencer(Sequencer *sequencerPtr)
{ sequencer=sequencerPtr; }
// start simulation over this Patch of atoms
void HomePatch::runSequencer(int numberOfCycles)
{ sequencer->run(numberOfCycles); }

void
HomePatch::readPatchMap() {
  PatchMap *p = PatchMap::Object();
  PatchID nnPatchID[PatchMap::MaxOneAway];

  patchMigrationCounter = numNeighbors 
    = PatchMap::Object()->oneAwayNeighbors(patchID, nnPatchID);
  DebugM( 1, "NumNeighbors for pid " <<patchID<<" is "<< numNeighbors << "\n");
  for (int n=0; n<numNeighbors; n++) {
    realInfo[n].destNodeID = p->node(realInfo[n].destPatchID = nnPatchID[n]);
     DebugM( 1, " nnPatchID=" <<nnPatchID[n]<<" nnNodeID="<< realInfo[n].destNodeID<< "\n");
    realInfo[n].mList = NULL;
  }

  // Make mapping from the 3x3x3 cube of pointers to real migration info
  for (int i=0; i<3; i++)
    for (int j=0; j<3; j++)
      for (int k=0; k<3; k++)
      {
	int pid =  p->pid(p->xIndex(patchID)+i-1, 
	    p->yIndex(patchID)+j-1, p->zIndex(patchID)+k-1);
	if (pid < 0) {
	   DebugM(5, "ERROR, for patchID " << patchID <<" I got neigh pid = " << pid << "\n");
	}
	if (pid == patchID) {
	  mInfo[i][j][k] = NULL;
	}
	else {
	  for (n = 0; n<numNeighbors; n++) {
	    if (pid == realInfo[n].destPatchID) {
	      mInfo[i][j][k] = &realInfo[n];
	      break;
	    }
	  }
	  if (n == numNeighbors) { // disaster! 
	    DebugM(4,"BAD News, I could not find PID " << pid << "\n");
	  }
	}
      }

  DebugM(1,"Patch("<<patchID<<") # of neighbors = " << numNeighbors << "\n");
}

HomePatch::~HomePatch()
{
}


void HomePatch::boxClosed(int)
{
  if ( ! --boxesOpen )
  {
    DebugM(1,patchID << ": " << CthSelf() << " awakening sequencer "
	<< sequencer->thread << "(" << patchID << ") @" << CmiTimer() << "\n");
    // only awaken suspended threads.  Then say it is suspended.
    sequencer->awaken();
    return;
  }
  else
  {
    DebugM(1,patchID << ": " << boxesOpen << " boxes left to close.\n");
  }
}

void HomePatch::registerProxy(RegisterProxyMsg *msg) {
  DebugM(4, "registerProxy("<<patchID<<") - adding node " <<msg->node<<"\n");
  proxy.add(ProxyListElem(msg->node,forceBox.checkOut()));
  ProxyAtomsMsg *nmsg = new (MsgIndex(ProxyAtomsMsg)) ProxyAtomsMsg;
  nmsg->patch = patchID;
  nmsg->atomIDList = atomIDList;
  nmsg->prepack();
  ProxyMgr::Object()->sendProxyAtoms(nmsg,msg->node);
  delete msg;
}

void HomePatch::unregisterProxy(UnregisterProxyMsg *msg) {
  int i = proxy.findIndex(ProxyListElem(msg->node));
  forceBox.checkIn(proxy[i].forceBox);
  proxy.del(i);
  delete msg;
}

void HomePatch::receiveResults(ProxyResultMsg *msg)
{
  DebugM(4, "patchID("<<patchID<<") receiveRes() nodeID("<<msg->node<<")\n");
  int i = proxy.findIndex(ProxyListElem(msg->node));
  Results *r = proxy[i].forceBox->open();
  for ( int k = 0; k < Results::maxNumForces; ++k )
  {
    Force *f = r->f[k];
    for ( register int j = 0; j < numAtoms; ++j )
    {
      f[j] += msg->forceList[k][j];
    }
  }
  proxy[i].forceBox->close(&r);
  delete msg;
}


void HomePatch::positionsReady(int doMigration)
{
  if (!patchMapRead) {
    readPatchMap();
    patchMapRead = 1;
  }
      
  doMigration = doMigration && numNeighbors;

  if (doMigration) {
    doAtomMigration();
  }

  // Must Add Proxy Changes when migration completed!
  ProxyListIter pli(proxy);
  for ( pli = pli.begin(); pli != pli.end(); ++pli )
  {
    if (doMigration) {
      ProxyAllMsg *allmsg 
	= new (MsgIndex(ProxyAllMsg)) ProxyAllMsg;
      allmsg->patch = patchID;
      allmsg->flags = flags;
      allmsg->positionList = p;
      allmsg->atomIDList = atomIDList;
      DebugM(1, "atomIDList.size() = " << atomIDList.size() << " p.size() = " << p.size() << "\n" );
      ProxyMgr::Object()->sendProxyAll(allmsg,pli->node);
    } else {
      ProxyDataMsg *nmsg 
	= new (MsgIndex(ProxyDataMsg)) ProxyDataMsg;
      nmsg->patch = patchID;
      nmsg->flags = flags;
      nmsg->positionList = p;
      ProxyMgr::Object()->sendProxyData(nmsg,pli->node);
    }   
  }
  DebugM(4, "patchID("<<patchID<<") doing positions Ready\n");
  Patch::positionsReady(doMigration);

}


void HomePatch::addForceToMomentum(const BigReal timestep, const int ftag)
{
  const BigReal dt = timestep / TIMEFACTOR;
  //if (v.check() == NULL || f.check() == NULL || a.check() == NULL) {
      //DebugM(5, "NULL found! v="<<v<<" f="<<f<<" a="<<a<<"\n");
  //}
  for ( int i = 0; i < numAtoms; ++i )
  {
    v[i] += f[ftag][i] * ( dt / a[i].mass );
    if ( a[i].flags & ATOM_FIXED ) v[i] = 0;
  }
}

void HomePatch::addVelocityToPosition(const BigReal timestep)
{
  const BigReal dt = timestep / TIMEFACTOR;
  for ( int i = 0; i < numAtoms; ++i )
  {
    if ( ! ( a[i].flags & ATOM_FIXED ) ) p[i] += v[i] * dt;
  }
}

BigReal HomePatch::calcKineticEnergy()
{
  BigReal total = 0;
  for ( int i = 0; i < numAtoms; ++i )
  {
     total += 0.5 * a[i].mass * v[i] * v[i];
  }
  return total;
}

Vector HomePatch::calcMomentum()
{
  Vector total;
  for ( int i = 0; i < numAtoms; ++i )
  {
     total += a[i].mass * v[i];
  }
  return total;
}

Vector HomePatch::calcAngularMomentum()
{
  Vector total;
  for ( int i = 0; i < numAtoms; ++i )
  {
     total += cross(a[i].mass,p[i],v[i]); // m r % v
  }
  return total;
}

void HomePatch::submitLoadStats(int timestep)
{
  LdbCoordinator::Object()->patchLoad(patchID,numAtoms,timestep);
}


void
HomePatch::doAtomMigration()
{
  int i;
  int xdev, ydev, zdev;
  MigrationList *mCur;

  // Drain the migration message buffer
  //for (i=0; i<numMlBuf; i++) {
  //   DebugM(3, "Draining migration buffer ("<<i<<","<<numMlBuf<<")\n");
  //   depositMigration(srcID[i], mlBuf[i]);
  //}
  //numMlBuf = 0;
     
  // realInfo points to migration lists for neighbors we actually have. 
  //    element of mInfo[3][3][3] points to an element of realInfo
  for (i=0; i<numNeighbors; i++) {
    realInfo[i].mList = NULL;
  }

  // Purge the AtomMap
  AtomMap::Object()->unregisterIDs(patchID,atomIDList);

  // Determine atoms that need to migrate
  i = 0;
  SimParameters *simParams = Node::Object()->simParameters;
  Molecule *mol = Node::Object()->molecule;
  while ( i < atomIDList.size() )
  {
     Position Min = lattice.unscale(min);
     Position Max = lattice.unscale(max);

/*
     if (simParams->splitPatch == SPLIT_PATCH_HYDROGEN)
	{
	// All atoms are always in hydrogen-group order.  Thus, if it's
	// a group member, then it is always moved *after* the group parent.
	// So, we only determine where to move if it's a group parent.

	// group members have their group parent ID
	// group parents have their own ID
	// ...do this in stages so it's easier to read
	j = mol->is_hydrogenGroupParent(atomIDList[i]);	// j is parent flag
	}
      else j=1;	// check all atoms
*/

      if ( a[i].hydrogenGroupSize )
	  {
	  // check if atom should is within bounds
	  if (p[i].x < Min.x) xdev = 0;
	  else if (Max.x <= p[i].x) xdev = 2; 
	  else xdev = 1;

	  if (p[i].y < Min.y) ydev = 0;
	  else if (Max.y <= p[i].y) ydev = 2; 
	  else ydev = 1;

	  if (p[i].z < Min.z) zdev = 0;
	  else if (Max.z <= p[i].z) zdev = 2; 
	  else zdev = 1;
	  }

     if (mInfo[xdev][ydev][zdev]) { // process atom for migration
                                    // Don't migrate if destination is myself

       // See if we have a migration list already
       if (NULL == (mCur = mInfo[xdev][ydev][zdev]->mList)) {
	 // new: all mList pointers are actually in realInfo[].mList
	 //one of the following below sendMigrationMsg() does delete
	 // [1] in pack of MigrateAtomsMsg (PatchMgr.C)
	 // [2] in recvMigrateAtoms (PatchMgr.C)
	 mCur = mInfo[xdev][ydev][zdev]->mList = new MigrationList;
       }
       Force force[Results::maxNumForces];
       for ( int j = 0; j < Results::maxNumForces; ++j ) force[j] = f[j][i];
       DebugM(3,"Migrating atom " << atomIDList[i] << " from patch "
		<< patchID << " with position " << p[i] << "\n");
       mCur->add(MigrationElem(atomIDList[i], a[i], p[i], p[i], v[i], force)
       );
       a.del(i);
       atomIDList.del(i,1);
       p.del(i);
       // pInit.del(i);
       v.del(i);
       for ( j = 0; j < Results::maxNumForces; ++j ) f[j].del(i);
       // f_short.del(i);
       // f_long.del(i);
     }
     else
     {
       ++i;
     }
  }
  numAtoms = atomIDList.size();

  PatchMgr::Object()->sendMigrationMsgs(patchID, realInfo, numNeighbors);
  // for (i=0; i < numNeighbors; i++) {
  //   PatchMgr::Object()->sendMigrationMsg(patchID, realInfo[i]);
  // }

  // signal depositMigration() that we are inMigration mode
  inMigration = true;

  // Drain the migration message buffer
  for (i=0; i<numMlBuf; i++) {
     DebugM(1, "Draining migration buffer ("<<i<<","<<numMlBuf<<")\n");
     depositMigration(msgbuf[i]);
  }
  numMlBuf = 0;
     
  if (!allMigrationIn) {
    DebugM(3,"All Migrations NOT in, we are suspending patch "<<patchID<<"\n");
    migrationSuspended = true;
    sequencer->suspend();
    migrationSuspended = false;
  }
  allMigrationIn = false;
  // indexAtoms();  NEVER USED -JCP

  // reload the AtomMap
  AtomMap::Object()->registerIDs(patchID,atomIDList);

  inMigration = false;
}

void 
HomePatch::depositMigration(MigrateAtomsMsg *msg)
{
  PatchID srcPatchID;
  MigrationList *migrationList;

  if (!inMigration) { // We have to buffer changes due to migration
		      // until our patch is in migration mode
    msgbuf[numMlBuf++] = msg;
    return;
  } 

  srcPatchID = msg->srcPatchID;
  migrationList = msg->migrationList;

  if (migrationList) {
    MigrationListIter mi(*migrationList);
    for (mi = mi.begin(); mi != mi.end(); mi++) {
      DebugM(1,"Migrating atom " << mi->atomID << " to patch "
		<< patchID << " with position " << mi->pos << "\n"); 
      a.add(mi->atomProp);
      atomIDList.add(mi->atomID);
      p.add(lattice.nearest(mi->pos,center));
      v.add(mi->vel);
      for ( int j = 0; j < Results::maxNumForces; ++j )
        f[j].add(mi-> force[j]);
    }
    delete migrationList;
    migrationList = NULL;
  }
  numAtoms = atomIDList.size();
  delete msg;

  DebugM(3,"Counter on " << patchID << " = " << patchMigrationCounter << "\n");
  if (!--patchMigrationCounter) {
    DebugM(3,"All Migrations are in for patch "<<patchID<<"\n");
    allMigrationIn = true;
    patchMigrationCounter = numNeighbors;
    if (migrationSuspended) {
      DebugM(3,"patch "<<patchID<<" is being awakened\n");
      migrationSuspended = false;
      sequencer->awaken();
      return;
    }
    else {
       DebugM(3,"patch "<<patchID<<" was not suspended\n");
    }
  }
}


/***************************************************************************
 * RCS INFORMATION:
 *
 *	$RCSfile: HomePatch.C,v $
 *	$Author: jim $	$Locker:  $		$State: Exp $
 *	$Revision: 1.1037 $	$Date: 1997/10/06 00:12:31 $
 *
 ***************************************************************************
 * REVISION HISTORY:
 *
 * $Log: HomePatch.C,v $
 * Revision 1.1037  1997/10/06 00:12:31  jim
 * Added PatchMap.inl, sped up cycle-boundary tuple code.
 *
 * Revision 1.1036  1997/09/28 10:19:07  milind
 * Fixed priorities, ReductionMgr etc.
 *
 * Revision 1.1035  1997/09/19 08:55:31  jim
 * Added rudimentary but relatively efficient fixed atoms.  New options
 * are fixedatoms, fixedatomsfile, and fixedatomscol (nonzero means fixed).
 * Energies will be affected, although this can be fixed with a little work.
 *
 * Revision 1.1034  1997/09/19 05:17:43  jim
 * Cleaned up and tweaked hydrogen-group based temporary pairlist
 * generation for roughly a 6% performance improvement.
 *
 * Revision 1.1033  1997/08/26 16:26:15  jim
 * Revamped prioritites for petter performance and easier changes.
 *
 * Revision 1.1032  1997/08/22 20:12:03  milind
 * Turned on Priorities.
 *
 * Revision 1.1031  1997/07/08 15:48:08  milind
 * Made namd2 to work with Origin2000: Again...
 *
 * Revision 1.1030  1997/04/21 00:58:33  jim
 * Fixed hang on patch migration in systems with only one patch.
 *
 * Revision 1.1029  1997/04/10 22:29:11  jim
 * First steps towards combining atom migration messages.
 *
 * Revision 1.1028  1997/04/10 09:13:57  ari
 * Final debugging for compute migration / proxy creation for load balancing.
 * Lots of debug code added, mostly turned off now.
 * Fixed bug in PositionBox when Patch had no dependencies.
 * Eliminated use of cout and misuse of iout in numerous places.
 *                                            Ari & Jim
 *
 * Revision 1.1027  1997/04/08 07:08:35  ari
 * Modification for dynamic loadbalancing - moving computes
 * Still bug in new computes or usage of proxies/homepatches.
 * Works if ldbStrategy is none as before.
 *
 * Revision 1.1026  1997/04/06 22:45:04  ari
 * Add priorities to messages.  Mods to help proxies without computes.
 * Added quick enhancement to end of list insertion of ResizeArray(s)
 *
 * Revision 1.1025  1997/03/31 16:12:51  nealk
 * Atoms now can migrate by hydrogen groups.
 *
 * Revision 1.1024  1997/03/27 20:25:44  brunner
 * Changes for LdbCoordinator, the load balance control BOC
 *
 * Revision 1.1023  1997/03/27 08:04:16  jim
 * Reworked Lattice to keep center of cell fixed during rescaling.
 *
 * Revision 1.1022  1997/03/18 18:09:01  jim
 * Revamped collection system to ensure ordering and eliminate
 * unnecessary collections.  Also reduced make dependencies.
 *
 * Revision 1.1021  1997/03/12 22:06:40  jim
 * First step towards multiple force returns and multiple time stepping.
 *
 * Revision 1.1020  1997/03/10 17:40:11  ari
 * UniqueSet changes - some more commenting and cleanup
 *
 * Revision 1.1019  1997/03/06 22:06:02  ari
 * Removed Compute.ci
 * Comments added - more code cleaning
 *
 * Revision 1.1018  1997/02/28 04:47:08  jim
 * Full electrostatics now works with fulldirect on one node.
 *
 * Revision 1.1017  1997/02/26 21:39:27  jim
 * Fixed migration with periodic boundary conditions to correctly
 * re-center migrated atoms on their new home patch.
 *
 * Revision 1.1016  1997/02/26 16:53:09  ari
 * Cleaning and debuging for memory leaks.
 * Adding comments.
 * Removed some dead code due to use of Quiescense detection.
 *
 * Revision 1.1015  1997/02/17 23:46:59  ari
 * Added files for cleaning up atom migration code
 *
 * Revision 1.1014  1997/02/13 23:17:17  ari
 * Fixed a final bug in AtomMigration - numatoms in ComputePatchPair.C not
 * set correctly in atomUpdate()
 *
 * Revision 1.1013  1997/02/13 17:06:22  jim
 * Turned off debugging.
 *
 * Revision 1.1012  1997/02/13 16:17:13  ari
 * Intermediate debuging commit - working to fix deep bug in migration?
 *
 * Revision 1.1011  1997/02/13 04:43:09  jim
 * Fixed initial hanging (bug in PatchMap, but it still shouldn't have
 * happened) and saved migration messages in the buffer from being
 * deleted, but migration still dies (even on one node).
 *
 * Revision 1.1010  1997/02/11 18:51:46  ari
 * Modified with #ifdef DPMTA to safely eliminate DPMTA codes
 * fixed non-buffering of migration msgs
 * Migration works on multiple processors
 *
 * Revision 1.1009  1997/02/11 16:31:48  nealk
 * Using a flag to determine suspend/awaken in sequencer.
 * Migration turned off (for now).
 *
 * Revision 1.1008  1997/02/10 19:44:29  nealk
 * More debugging code.  Corrected HomePatch/Sequencer awaken/suspend bug.
 *
 * Revision 1.1007  1997/02/10 08:17:30  jim
 * Commented out sending and allocation of unused data.
 *
 * Revision 1.1006  1997/02/07 17:39:38  ari
 * More debugging for atomMigration.
 * Using -w on CC got us some minor fixes
 * using purify got us a major memory problem due to bad sizing of dummy force
 *
 * Revision 1.1005  1997/02/07 07:51:42  jim
 * pInit aliases p, leading to problems.  I commented out any add or del
 * involving it.  Now seems to integrate correctly.
 *
 * Revision 1.1004  1997/02/07 05:42:30  ari
 * Some bug fixing - atom migration on one node works
 * Atom migration on multiple nodes gets SIGSEGV
 *
 * Revision 1.1003  1997/02/06 23:25:07  jim
 * Fixed bugs.
 *
 * Revision 1.1002  1997/02/06 21:20:50  jim
 * Fixed a couple of atom migration bugs.
 *
 * Revision 1.1001  1997/02/06 18:05:28  nealk
 * Modified (added some, turned off others) debug statements.
 *
 * Revision 1.1000  1997/02/06 15:58:26  ari
 * Resetting CVS to merge branches back into the main trunk.
 * We will stick to main trunk development as suggested by CVS manual.
 * We will set up tags to track fixed points of development/release
 * as suggested by CVS manual - all praise the CVS manual.
 *
 * Revision 1.779  1997/02/06 15:53:11  ari
 * Updating Revision Line, getting rid of branches
 *
 * Revision 1.778.2.1  1997/02/05 22:18:13  ari
 * Added migration code - Currently the framework is
 * there with compiling code.  This version does
 * crash shortly after migration is complete.
 * Migration appears to complete, but Patches do
 * not appear to be left in a correct state.
 *
 * Revision 1.778  1997/01/28 00:30:35  ari
 * internal release uplevel to 1.778
 *
 * Revision 1.777.2.2  1997/01/27 22:45:14  ari
 * Basic Atom Migration Code added.
 * Added correct magic first line to .h files for xemacs to go to C++ mode.
 * Compiles and runs without migration turned on.
 *
 * Revision 1.777.2.1  1997/01/21 23:04:44  ari
 * Basic framework for atom migration placed into code.  - Non
 * functional since it is not called.  Works currently without
 * atom migration.
 *
 * Revision 1.777  1997/01/17 19:36:10  ari
 * Internal CVS leveling release.  Start development code work
 * at 1.777.1.1.
 *
 * Revision 1.19  1997/01/15 17:59:11  jim
 * now divide by TIMEFACTOR
 *
 * Revision 1.18  1997/01/15 17:09:43  ari
 * Commented out new Atom migration code
 *
 * Revision 1.17  1997/01/10 22:38:37  jim
 * kinetic energy reporting
 *
 * Revision 1.16  1996/12/17 23:58:02  jim
 * proxy result reporting is working
 *
 * Revision 1.15  1996/12/17 22:13:22  jim
 * implemented ProxyDataMsg use
 *
 * Revision 1.14  1996/12/17 08:55:25  jim
 * added node argument to sendProxyAtoms
 *
 * Revision 1.13  1996/12/16 22:52:43  jim
 * added placement new and explicit destructor calls to ProxyAtomsMsg
 *
 * Revision 1.12  1996/12/14 00:02:42  jim
 * debugging ProxyAtomsMsg path to make compute creation work
 *
 * Revision 1.11  1996/12/11 22:31:41  jim
 * added integration methods for Sequencer
 *
 * Revision 1.10  1996/12/05 23:45:09  ari
 * *** empty log message ***
 *
 * Revision 1.9  1996/12/05 01:44:16  ari
 * started toward proxy management
 *
 * Revision 1.8  1996/12/01 02:35:47  jim
 * added patchID reporting
 *
 * Revision 1.7  1996/11/30 00:35:51  jim
 * implemented boxClosed(), useSequencer(), runSequencer()
 *
 * Revision 1.6  1996/11/22 01:44:53  jim
 * added calls to service AtomMap
 *
 * Revision 1.5  1996/10/16 08:22:39  ari
 * *** empty log message ***
 *
 * Revision 1.4  1996/10/04 21:14:47  jim
 * Moved functionality to Patch
 *
 * Revision 1.3  1996/09/03 22:54:09  ari
 * *** empty log message ***
 *
 * Revision 1.2  1996/08/29 00:50:42  ari
 * *** empty log message ***
 *
 * Revision 1.1  1996/08/19 22:07:49  ari
 * Initial revision
 *
 * Revision 1.7  1996/07/15 21:18:15  gursoy
 * *** empty log message ***
 *
 * Revision 1.6  1996/07/10 16:25:29  gursoy
 * *** empty log message ***
 *
 * Revision 1.5  1996/07/10 16:19:52  gursoy
 * *** empty log message ***
 *
 * Revision 1.4  1996/07/02 15:10:27  gursoy
 * *** empty log message ***
 *
 * Revision 1.3  1996/06/11 22:36:23  brunner
 * *** empty log message ***
 *
 * Revision 1.2  1996/06/11 20:07:22  gursoy
 * *** empty log message ***
 *
 * Revision 1.1  1996/05/30 21:31:36  gursoy
 * Initial revision
 *
 ***************************************************************************/
