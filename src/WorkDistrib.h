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

const double patchSize = 4.0;

enum { maxPatchDepends = 126 };

class MapDistribMsg;

class WorkDistrib : public BOCclass
{
public:
  WorkDistrib(InitMsg *msg);
  ~WorkDistrib(void);

  static void messageMovePatchDone();
  void movePatchDone(DoneMsg *msg);

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
    char *computeMapData = (char*)computeMap->pack(&computeMapSize);
    *length = sizeof(int) + patchMapSize + sizeof(int) + computeMapSize;
    char *buffer = (char*)new_packbuffer(this,*length);
    char *b = buffer;
    *((int*)b) = patchMapSize;
    b += sizeof(int);
    memcpy(b,patchMapData,patchMapSize);
    delete [] patchMapData;
    b += computeMapSize;
    *((int*)b) = computeMapSize;
    b += sizeof(int);
    memcpy(b,computeMapData,computeMapSize);
    delete [] computeMapData;
    return buffer;
  }
  void unpack (void *in)
  {
    char *buffer = (char*)in;
    int patchMapSize = *((int*)buffer);
    buffer += sizeof(int);
    patchMap = PatchMap::Instance();
    patchMap->unpack((void*)buffer);
    buffer += patchMapSize;
    int computeMapSize = *((int*)buffer);
    buffer += sizeof(int);
    computeMap = ComputeMap::Instance();
    computeMap->unpack((void*)buffer);
  }
};


#endif /* WORKDISTRIB_H */

/***************************************************************************
 * RCS INFORMATION:
 *
 *	$RCSfile: WorkDistrib.h,v $
 *	$Author: jim $	$Locker:  $		$State: Exp $
 *	$Revision: 1.14 $	$Date: 1996/12/12 08:57:17 $
 *
 ***************************************************************************
 * REVISION HISTORY:
 *
 * $Log: WorkDistrib.h,v $
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
