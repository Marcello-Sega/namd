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

static char ident[] = "@(#)$Header: /home/cvs/namd/cvsroot/namd2/src/HomePatch.C,v 1.1006 1997/02/07 17:39:38 ari Exp $";

#include "ckdefs.h"
#include "chare.h"
#include "c++interface.h"

#include "HomePatch.h"
#include "AtomMap.h"
#include "PatchMap.h"
#include "main.h"
#include "ProxyMgr.top.h"
#include "ProxyMgr.h"
#include "Migration.h"
#include "PatchMgr.h"

#define MIN_DEBUG_LEVEL 3
// #define DEBUGM
#include "Debug.h"

HomePatch::HomePatch(PatchID pd, AtomIDList al, PositionList pl, 
		     VelocityList vl) : Patch(pd,al,pl), pInit(pl), v(vl) 
{ 
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

  migrationSuspended = false;
  allMigrationIn = false;
  patchMapRead = 0; // We delay read of PatchMap data
		    // to make sure it is really valid
}

void
HomePatch::readPatchMap() {
  PatchMap *p = PatchMap::Object();
  PatchID nnPatchID[PatchMap::MaxOneAway];

  patchMigrationCounter = numNeighbors 
    = PatchMap::Object()->oneAwayNeighbors(patchID, nnPatchID);
  DebugM( 4, "NumNeighbors for pid " <<patchID<<" is "<< numNeighbors << "\n");
  for (int n=0; n<numNeighbors; n++) {
    realInfo[n].destNodeID = p->node(realInfo[n].destPatchID = nnPatchID[n]);
     DebugM( 4, " nnPatchID=" <<nnPatchID[n]<<" nnNodeID="<< realInfo[n].destNodeID<< "\n");
    realInfo[n].mList = NULL;
  }

  // Make mapping from the 3x3x3 cube of pointers to real migration info
  for (int i=0; i<3; i++)
    for (int j=0; j<3; j++)
      for (int k=0; k<3; k++)
      {
	int pid =  p->pid(p->xIndex(patchID)+i-1, 
	    p->yIndex(patchID)+j-1, p->zIndex(patchID)+k-1);
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

  DebugM(4,"Patch("<<patchID<<") # of neighbors = " << numNeighbors << "\n");
}

HomePatch::~HomePatch()
{
}


void HomePatch::boxClosed(int)
{
  if ( ! --boxesOpen )
  {
    DebugM(2,patchID << ": " << "Trying to awaken sequencer.\n");
    sequencer->awaken();
  }
  else
  {
    DebugM(2,patchID << ": " << boxesOpen << " boxes left to close.\n");
  }
}

void HomePatch::registerProxy(RegisterProxyMsg *msg) {
  proxy.add(ProxyListElem(msg->node,forceBox.checkOut()));
  ProxyAtomsMsg *nmsg = new (MsgIndex(ProxyAtomsMsg)) ProxyAtomsMsg;
  nmsg->patch = patchID;
  nmsg->atomIDList = atomIDList;
  ProxyMgr::Object()->sendProxyAtoms(nmsg,msg->node);
}

void HomePatch::unregisterProxy(UnregisterProxyMsg *msg) {
  int i = proxy.findIndex(ProxyListElem(msg->node));
  forceBox.checkIn(proxy[i].forceBox);
  proxy.del(i);
}

void HomePatch::receiveResults(ProxyResultMsg *msg)
{
  int i = proxy.findIndex(ProxyListElem(msg->node));
  Force* f = proxy[i].forceBox->open();
  for ( int j = 0; j < numAtoms; ++j )
  {
    f[j] += msg->forceList[j];
  }
  proxy[i].forceBox->close(&f);
}


void HomePatch::positionsReady(int doMigration)
{
  if (!patchMapRead) {
    readPatchMap();
    patchMapRead = 1;
  }
  if (doMigration) {
    doAtomMigration();
  }

  // Must Add Proxy Changes when migration completed!
  ProxyListIter pli(proxy);
  for ( pli = pli.begin(); pli != pli.end(); ++pli )
  {
    if (doMigration) {
      ProxyAllMsg *allmsg = new (MsgIndex(ProxyAllMsg)) ProxyAllMsg;
      allmsg->patch = patchID;
      allmsg->positionList = p;
      allmsg->atomIDList = atomIDList;
      DebugM(4, "atomIDList.size() = " << atomIDList.size() 
	<< " positionList.size() = " << positionList.size() << "\n" );
      ProxyMgr::Object()->sendProxyAll(allmsg,pli->node);
    } else {
      ProxyDataMsg *nmsg = new (MsgIndex(ProxyDataMsg)) ProxyDataMsg;
      nmsg->patch = patchID;
      nmsg->positionList = p;
      ProxyMgr::Object()->sendProxyData(nmsg,pli->node);
    }   
  }
  Patch::positionsReady(doMigration);
}


void HomePatch::addForceToMomentum(const BigReal timestep)
{
  const BigReal dt = timestep / TIMEFACTOR;
  for ( int i = 0; i < numAtoms; ++i )
  {
    v[i] += f[i] * ( dt / a[i].mass );
  }
}

void HomePatch::addVelocityToPosition(const BigReal timestep)
{
  const BigReal dt = timestep / TIMEFACTOR;
  for ( int i = 0; i < numAtoms; ++i )
  {
    p[i] += v[i] * dt;
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


void
HomePatch::doAtomMigration()
{
  int i;
  int xdev, ydev, zdev;
  MigrationList *mCur;

  // Null pointers to migration lists.
  for (i=0; i<numNeighbors; i++) {
    realInfo[i].mList = NULL;
  }

  // Purge the AtomMap
  AtomMap::Object()->unregisterIDs(patchID,atomIDList);

  // Determine atoms that need to migrate
  i = 0;
  while ( i < atomIDList.size() )
  {
     if (p[i].x < min.x) xdev = 0;
     else if (max.x <= p[i].x) xdev = 2; 
     else xdev = 1;

     if (p[i].y < min.y) ydev = 0;
     else if (max.y <= p[i].y) ydev = 2; 
     else ydev = 1;

     if (p[i].z < min.z) zdev = 0;
     else if (max.z <= p[i].z) zdev = 2; 
     else zdev = 1;

     if (mInfo[xdev][ydev][zdev]) { // process atom for migration
                                    // Don't migrate if destination is myself

       // See if we have a migration list already
       if (NULL == (mCur = mInfo[xdev][ydev][zdev]->mList)) {
	 mCur = mInfo[xdev][ydev][zdev]->mList = new MigrationList;
       }
       DebugM(4,"Migrating atom " << atomIDList[i] << " from patch "
		<< patchID << " with position " << p[i] << "\n");
       mCur->add(MigrationElem(atomIDList[i], a[i], pInit[i],
         p[i], v[i], f[i], f_short[i], f_long[i])
       );
       a.del(i);
       atomIDList.del(i,1);
       p.del(i);
       // pInit.del(i);
       v.del(i);
       f.del(i);
       f_short.del(i);
       f_long.del(i);
     }
     else
     {
       ++i;
     }
  }
  numAtoms = atomIDList.size();

  for (i=0; i < numNeighbors; i++) {
    PatchMgr::Object()->sendMigrationMsg(patchID, realInfo[i]);
  }

  if (!allMigrationIn) {
    DebugM(4,"All Migrations NOT in, we are suspending patch "<<patchID<<"\n");
    migrationSuspended = true;
    sequencer->suspend();
    migrationSuspended = false;
  }
  allMigrationIn = false;
  indexAtoms();

  // reload the AtomMap
  AtomMap::Object()->registerIDs(patchID,atomIDList);
}

void 
HomePatch::depositMigration(PatchID srcPatchID, MigrationList *migrationList)
{
  DebugM(4,"depositMigration from "<<srcPatchID<<" on "<<patchID<<"\n");
  if (migrationList) {
    MigrationListIter mi(*migrationList);
    for (mi = mi.begin(); mi != mi.end(); mi++) {
      DebugM(4,"Migrating atom " << mi->atomID << " to patch "
		<< patchID << " with position " << mi->pos << "\n"); 
      a.add(mi->atomProp);
      atomIDList.add(mi->atomID);
      p.add(mi->pos);
      // pInit.add(mi->posInit);
      v.add(mi->vel);
      f.add(mi->force);
      f_short.add(mi->forceShort);
      f_long.add(mi->forceLong);
    }
  }
  numAtoms = atomIDList.size();

  DebugM(4,"Counter on " << patchID << " = " << patchMigrationCounter << "\n");
  if (!--patchMigrationCounter) {
    DebugM(4,"All Migrations are in for patch "<<patchID<<"\n");
    allMigrationIn = true;
    patchMigrationCounter = numNeighbors;
    if (migrationSuspended) {
      DebugM(4,"patch "<<patchID<<" is being awakened\n");
      migrationSuspended = false;
      sequencer->awaken();
    }
  }
}


/***************************************************************************
 * RCS INFORMATION:
 *
 *	$RCSfile: HomePatch.C,v $
 *	$Author: ari $	$Locker:  $		$State: Exp $
 *	$Revision: 1.1006 $	$Date: 1997/02/07 17:39:38 $
 *
 ***************************************************************************
 * REVISION HISTORY:
 *
 * $Log: HomePatch.C,v $
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
