//-*-c++-*-
/***************************************************************************/
/*                                                                         */
/*              (C) Copyright 1996 The Board of Trustees of the            */
/*                          University of Illinois                         */
/*                           All Rights Reserved                           */
/*									   */
/***************************************************************************/

/***************************************************************************/
/* DESCRIPTION:								   */
/*									   */
/***************************************************************************/

#ifndef _WORKDISTRIB_H
#define _WORKDISTRIB_H

#include "charm++.h"

#include "main.h"

#include "NamdTypes.h"
#include "BOCgroup.h"
#include "ComputeMap.h"
#include "WorkDistrib.decl.h"

class Node;
class LocalWorkMsg;
class Molecule;

enum { maxPatchDepends = 126 };

class MapDistribMsg;
class ComputeMapChangeMsg;
class ComputeMapDistribMsg;

class WorkDistrib : public BOCclass
{
public:
  WorkDistrib();
  ~WorkDistrib(void);

  // static void messageMovePatchDone();
  // void movePatchDone();

  static void messageEnqueueWork(Compute *);
  void enqueueWork(LocalWorkMsg *msg);
  void enqueueSelfA(LocalWorkMsg *msg);
  void enqueueSelfB(LocalWorkMsg *msg);
  void enqueueWorkA(LocalWorkMsg *msg);
  void enqueueWorkB(LocalWorkMsg *msg);
  void enqueueWorkC(LocalWorkMsg *msg);

  void mapComputes(void);
  void sendMaps(void);
  void saveComputeMapChanges(int,int);
  void recvComputeMapChanges(ComputeMapChangeMsg *);
  void doneSaveComputeMap();
  void createHomePatches(void);
  void distributeHomePatches(void);
  void patchMapInit(void);
  void assignNodeToPatch(void);

  void saveMaps(MapDistribMsg *msg);

private:
  void mapComputeNonbonded(void);
  void mapComputeHomePatches(ComputeType);
  void mapComputePatch(ComputeType);
  void assignPatchesToLowestLoadNode(void);
  void assignPatchesRecursiveBisection(void);
  void assignPatchesRoundRobin(void);
  void velocities_from_PDB(char *filename, 
			   Vector *v, int totalAtoms);
  void velocities_from_binfile(char *fname, Vector *vels, int n);
  void random_velocities(BigReal Temp, Molecule *structure,
			 Vector *v, int totalAtoms);
  void remove_com_motion(Vector *vel, Molecule *structure, int n);

  Boolean mapsArrived;
  Boolean awaitingMaps;
  CthThread awaitingMapsTh;

  int saveComputeMapReturnEP;
  int saveComputeMapReturnChareID;
  int saveComputeMapCount;
};

#include <string.h>
#include "PatchMap.h"
#include "ComputeMap.h"

class ComputeMapDistribMsg : public CMessage_ComputeMapDistribMsg
{
public:
  ComputeMapDistribMsg(void) : computeMap(0) { ; }
  ComputeMap *computeMap;

  // pack and unpack functions
  static void* pack(ComputeMapDistribMsg* msg)
  {
    int computeMapSize = msg->computeMap->packSize();
    char *buffer = (char*)CkAllocBuffer(msg,computeMapSize);
    msg->computeMap->pack(buffer);
    delete msg;
    return buffer;
  }

  static ComputeMapDistribMsg* unpack(void *ptr)
  {
    void *_ptr = CkAllocBuffer(ptr, sizeof(ComputeMapDistribMsg));
    ComputeMapDistribMsg *m = new (_ptr) ComputeMapDistribMsg();
    m->computeMap = ComputeMap::Object();
    m->computeMap->unpack((char*)ptr);
    CkFreeMsg(ptr);
    return m;
  }
};

class MapDistribMsg : public CMessage_MapDistribMsg
{
public:
  MapDistribMsg(void) : patchMap(0), computeMap(0) { ; }
  PatchMap *patchMap;
  ComputeMap *computeMap;

  // pack and unpack functions
  static void* pack(MapDistribMsg *msg)
  {
    int patchMapSize = msg->patchMap->packSize();
    int computeMapSize = msg->computeMap->packSize();
    int length = sizeof(int) + patchMapSize + sizeof(int) + computeMapSize;
    char *buffer = (char*)CkAllocBuffer(msg,length);
    char *b = buffer;
    *((int*)b) = patchMapSize;
    b += sizeof(int);
    msg->patchMap->pack(b);
    b += patchMapSize;
    *((int*)b) = computeMapSize;
    b += sizeof(int);
    msg->computeMap->pack(b);
    delete msg;
    return buffer;
  }

  static MapDistribMsg* unpack(void *ptr)
  {
    void *_ptr = CkAllocBuffer(ptr, sizeof(MapDistribMsg));
    MapDistribMsg *m = new (_ptr) MapDistribMsg;
    char *buffer = (char*)ptr;
    int patchMapSize = *((int*)buffer);
    buffer += sizeof(int);
    m->patchMap = PatchMap::Object();
    if ( ! ( m->patchMap->patchData ) )
    {
      m->patchMap->unpack(buffer);
    }
    buffer += patchMapSize;
    // int computeMapSize = *((int*)buffer);
    buffer += sizeof(int);
    m->computeMap = ComputeMap::Object();
    if ( ! ( m->computeMap->computeData ) )
    {
      m->computeMap->unpack(buffer);
    }
    CkFreeMsg(ptr);
    return m;
  }
};

class ComputeMapChangeMsg : public CMessage_ComputeMapChangeMsg
{
public:
  int newNodes[20000];
  int numNewNodes;
};


#endif /* WORKDISTRIB_H */

/***************************************************************************
 * RCS INFORMATION:
 *
 *	$RCSfile: WorkDistrib.h,v $
 *	$Author: brunner $	$Locker:  $		$State: Exp $
 *	$Revision: 1.1016 $	$Date: 1999/05/11 23:56:54 $
 *
 ***************************************************************************
 * REVISION HISTORY:
 *
 * $Log: WorkDistrib.h,v $
 * Revision 1.1016  1999/05/11 23:56:54  brunner
 * Changes for new charm version
 *
 * Revision 1.1015  1998/10/24 19:58:05  jim
 * Eliminated warnings generated by g++ -Wall.
 *
 * Revision 1.1014  1998/07/03 20:09:56  brunner
 * Self-compute spliting creation changes.  I hope this works.
 *
 * Revision 1.1013  1998/07/02 21:00:03  brunner
 * Changed initial patch distribution, should work on more PES
 *
 * Revision 1.1012  1998/03/03 23:05:31  brunner
 * Changed include files for new simplified Charm++ include file structure.
 *
 * Revision 1.1011  1998/01/22 20:11:06  brunner
 * Modified the ComputeMap redistribution to send only new patch assignments.
 *
 * Revision 1.1010  1997/08/20 23:27:43  jim
 * Created multiple enqueueWork entry points to aid analysis.
 *
 * Revision 1.1009  1997/04/08 07:09:05  ari
 * Modification for dynamic loadbalancing - moving computes
 * Still bug in new computes or usage of proxies/homepatches.
 * Works if ldbStrategy is none as before.
 *
 * Revision 1.1008  1997/04/07 21:09:59  brunner
 * Added RecBisection for initial patch distrib
 *
 * Revision 1.1007  1997/03/19 05:50:16  jim
 * Added ComputeSphericalBC, cleaned up make dependencies.
 *
 * Revision 1.1006  1997/03/15 22:15:36  jim
 * Added ComputeCylindricalBC.  Doesn't break anything but untested and
 * cylinder is along x axis (will fix soon).
 *
 * Revision 1.1005  1997/03/14 21:40:17  ari
 * Reorganized startup to make possible inital load
 * balancing by changing methods in WorkDistrib.
 * Also made startup more transparent and easier
 * to modify.
 *
 * Revision 1.1004  1997/03/06 22:18:17  brunner
 * Made utility functions private functions of WorkDistrib.
 *
 * Revision 1.1003  1997/02/26 16:53:20  ari
 * Cleaning and debuging for memory leaks.
 * Adding comments.
 * Removed some dead code due to use of Quiescense detection.
 *
 * Revision 1.1002  1997/02/14 20:24:58  jim
 * Patches are now sized according to config file.
 *
 * Revision 1.1001  1997/02/14 19:07:32  nealk
 * Added new/delete comments.
 * Played with DPMTA.
 *
 * Revision 1.1000  1997/02/06 15:59:28  ari
 * Resetting CVS to merge branches back into the main trunk.
 * We will stick to main trunk development as suggested by CVS manual.
 * We will set up tags to track fixed points of development/release
 * as suggested by CVS manual - all praise the CVS manual.
 *
 * Revision 1.778  1997/01/28 00:31:31  ari
 * internal release uplevel to 1.778
 *
 * Revision 1.777.2.2  1997/01/27 22:45:44  ari
 * Basic Atom Migration Code added.
 * Added correct magic first line to .h files for xemacs to go to C++ mode.
 * Compiles and runs without migration turned on.
 *
 * Revision 1.777.2.1  1997/01/22 21:42:15  jim
 * Larger patches, no two-away computes, small tweak to inner loop.
 *
 * Revision 1.777  1997/01/17 19:37:06  ari
 * Internal CVS leveling release.  Start development code work
 * at 1.777.1.1.
 *
 * Revision 1.17  1996/12/14 00:01:59  jim
 * removed debug from .h file
 *
 * Revision 1.16  1996/12/13 22:58:01  nealk
 * Found pack/unpack bug and corrected it.  (wrong offset!)
 *
 * Revision 1.15  1996/12/13 08:53:49  jim
 * now moves patches
 *
 * Revision 1.14  1996/12/12 08:57:17  jim
 * added MapDistribMsg packing / unpacking routines
 *
 * Revision 1.13  1996/12/01 21:02:37  jim
 * now adds all existing compute objects to map
 *
 * Revision 1.12  1996/11/22 00:18:51  ari
 * *** empty log message ***
 *
 * Revision 1.11  1996/10/29 23:35:27  ari
 * *** empty log message ***
 *
 * Revision 1.10  1996/10/29 17:58:07  brunner
 * Did some stuff.  I forget what
 *
 * Revision 1.9  1996/10/16 08:22:39  ari
 * *** empty log message ***
 *
 * Revision 1.8  1996/10/04 22:23:50  brunner
 * Added createComputes
 *
 * Revision 1.7  1996/08/29 00:50:42  ari
 * *** empty log message ***
 *
 * Revision 1.6  1996/08/23 22:03:52  brunner
 * *** empty log message ***
 *
 * Revision 1.5  1996/08/19 21:37:02  brunner
 * Create Patches from PDB data
 *
 * Revision 1.4  1996/08/19 20:39:11  brunner
 * *** empty log message ***
 *
 * Revision 1.3  1996/08/16 21:41:11  brunner
 * *** empty log message ***
 *
 * Revision 1.2  1996/08/16 21:16:04  ari
 * *** empty log message ***
 *
 * Revision 1.1  1996/08/16 20:34:23  brunner
 * Initial revision
 *
 * Revision 1.6  1996/08/03 20:08:09  brunner
 * *** empty log message ***
 *
 * Revision 1.5  1996/07/01 16:25:30  brunner
 * *** empty log message ***
 *
 * Revision 1.4  1996/06/14 18:40:54  brunner
 * *** empty log message ***
 *
 * Revision 1.3  1996/06/11 22:38:04  brunner
 * *** empty log message ***
 *
 * Revision 1.2  1996/06/04 16:12:18  brunner
 * Create dummy patches
 *
 * Revision 1.1  1996/05/30 20:16:09  brunner
 * Initial revision
 *
 * Revision 1.1  1996/05/30 20:11:09  brunner
 * Initial revision
 *
 ***************************************************************************/
