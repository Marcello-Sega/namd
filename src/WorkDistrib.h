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

#include "ckdefs.h"
#include "chare.h"
#include "c++interface.h"
#include "main.h"

#include "NamdTypes.h"
#include "BOCgroup.h"

class Node;
class DoneMsg;
class LocalWorkMsg;
class Molecule;

enum { maxPatchDepends = 126 };

class MapDistribMsg;

class WorkDistrib : public BOCclass
{
public:
  WorkDistrib(InitMsg *msg);
  ~WorkDistrib(void);

  // static void messageMovePatchDone();
  // void movePatchDone(DoneMsg *msg);

  static void messageEnqueueWork(Compute *);
  void enqueueWork(LocalWorkMsg *msg); // This is for testing

  void buildMaps(void);
  void sendMaps(void);
  void createPatches(void);
  void createComputes(void);

  void saveMaps(MapDistribMsg *msg);
  void awaitMaps(void);

private:
  void mapPatches(void);
  void mapComputes(void);
  void mapComputeNonbonded(void);
  void mapComputeHomePatches(ComputeType);
  void velocities_from_PDB(char *filename, 
			   Vector *v, int totalAtoms);
  void velocities_from_binfile(char *fname, Vector *vels, int n);
  void random_velocities(BigReal Temp, Molecule *structure,
			 Vector *v, int totalAtoms);
  void remove_com_motion(Vector *vel, Molecule *structure, int n);

  Boolean mapsArrived;
  Boolean awaitingMaps;
  CthThread awaitingMapsTh;
};

#include <string.h>
#include "PatchMap.h"
#include "ComputeMap.h"

class MapDistribMsg : public comm_object
{
public:
  MapDistribMsg(void) : patchMap(0), computeMap(0) { ; }
  PatchMap *patchMap;
  ComputeMap *computeMap;

  // pack and unpack functions
  void * pack (int *length)
  {
    int patchMapSize, computeMapSize;
    char *patchMapData = (char*)patchMap->pack(&patchMapSize);
    char *computeMapData = (char*)computeMap->pack(&computeMapSize);	// "new" for data
    *length = sizeof(int) + patchMapSize + sizeof(int) + computeMapSize;
    char *buffer = (char*)new_packbuffer(this,*length);
    char *b = buffer;
    *((int*)b) = patchMapSize;
    b += sizeof(int);
    memcpy(b,patchMapData,patchMapSize);
    delete [] patchMapData;
    b += patchMapSize;
    *((int*)b) = computeMapSize;
    b += sizeof(int);
    memcpy(b,computeMapData,computeMapSize);
    delete [] computeMapData;	// allocated in ComputeMap.C::pack
    return buffer;
  }
  void unpack (void *in)
  {
    char *buffer = (char*)in;
    int patchMapSize = *((int*)buffer);
    buffer += sizeof(int);
    patchMap = PatchMap::Object();
    if ( ! ( patchMap->patchData ) )
    {
      patchMap->unpack((void*)buffer);
    }
    buffer += patchMapSize;
    int computeMapSize = *((int*)buffer);
    buffer += sizeof(int);
    computeMap = ComputeMap::Object();
    if ( ! ( computeMap->computeData ) )
    {
      computeMap->unpack((void*)buffer);
    }
  }
};


#endif /* WORKDISTRIB_H */

/***************************************************************************
 * RCS INFORMATION:
 *
 *	$RCSfile: WorkDistrib.h,v $
 *	$Author: brunner $	$Locker:  $		$State: Exp $
 *	$Revision: 1.1004 $	$Date: 1997/03/06 22:18:17 $
 *
 ***************************************************************************
 * REVISION HISTORY:
 *
 * $Log: WorkDistrib.h,v $
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
