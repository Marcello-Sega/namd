/***************************************************************************/
/*         (C) Copyright 1996,1997 The Board of Trustees of the            */
/*                          University of Illinois                         */
/*                           All Rights Reserved                           */
/***************************************************************************/
/***************************************************************************
 * DESCRIPTION: Methods are primarily for pack(ing) and unpack(ing)
 * 		messages for Charm.
 *		
 ***************************************************************************/

#include "charm++.h"

#include "Migration.h"
#include "MigrateAtomsMsg.h"
#include "NamdTypes.h"
// #define DEBUGM
#define MIN_DEBUG_LEVEL 3
#include "Debug.h"

#include "PatchMgr.decl.h"
#include "PatchMap.h"
#include "HomePatch.h"

MigrateAtomsMsg::MigrateAtomsMsg(void) { 
    migrationList = NULL; 
}

MigrateAtomsMsg::~MigrateAtomsMsg(void) { 
  // delete migrationList; NOW DELETED on pack() and by HomePatch on use!
}

MigrateAtomsMsg::MigrateAtomsMsg(PatchID src, PatchID dest, MigrationList *m) : 
      srcPatchID(src), destPatchID(dest), migrationList(m) 
{
    fromNodeID = CkMyPe();
}

void* MigrateAtomsMsg::pack (MigrateAtomsMsg *m) {
  int length = sizeof(NodeID) + sizeof(PatchID) + sizeof(PatchID) + sizeof(int) 
    + (m->migrationList != NULL 
	   ? m->migrationList->size() * sizeof(MigrationElem) 
	   : 0 );

  // This is what I want
  char *buffer;
  char *b = buffer = (char*)CkAllocBuffer(m,length);

  *((NodeID*)b) = m->fromNodeID; b += sizeof(NodeID);
  *((PatchID*)b) = m->srcPatchID; b += sizeof(PatchID);
  *((PatchID*)b) = m->destPatchID; b += sizeof(PatchID);

  if (m->migrationList != NULL) {
    *((int*)b) = m->migrationList->size(); b += sizeof(int);
    for ( int i = 0; i < m->migrationList->size(); i++ ) {
      *((MigrationElem*)b) = (*(m->migrationList))[i]; 
      b += sizeof(MigrationElem);
    }
  }
  else {
    *((int*)b) = 0; b += sizeof(int);
  }

  // Delete of new'd item from HomePatch::doAtomMigration()
  delete m->migrationList;
  m->migrationList = 0;

  delete m;
  return buffer;
}

MigrateAtomsMsg* MigrateAtomsMsg::unpack(void *ptr) {
  void *_ptr = CkAllocBuffer(ptr, sizeof(MigrateAtomsMsg));
  MigrateAtomsMsg* m = new (_ptr) MigrateAtomsMsg;
  char *b = (char*)ptr;
  m->fromNodeID = *((NodeID*)b); b += sizeof(NodeID);
  m->srcPatchID = *((PatchID*)b); b += sizeof(PatchID);
  m->destPatchID = *((PatchID*)b); b += sizeof(PatchID);
  int size = *((int*)b); b += sizeof(int);
  if (size != 0) {
    m->migrationList = new MigrationList();
    m->migrationList->resize(size);
    for ( int i = 0; i < size; i++ ) {
      (*(m->migrationList))[i] = *((MigrationElem*)b); 
      b += sizeof(MigrationElem);
    }
  }
  else {
    m->migrationList = NULL;
  }
  CkFreeMsg(ptr);
  return m;
}

MigrateAtomsCombinedMsg::MigrateAtomsCombinedMsg(void)
{
  fromNodeID = CkMyPe();
  totalAtoms = 0;
}

void MigrateAtomsCombinedMsg::
	add(PatchID source, PatchID destination, MigrationList *m)
{
  srcPatchID.add(source);
  destPatchID.add(destination);
  int n = ( m ? m->size() : 0 );
  numAtoms.add(n);
  totalAtoms += n;
  for ( int i = 0; i < n; ++i )
  {
    migrationList.add((*m)[i]);
  }
  delete m;
}


void MigrateAtomsCombinedMsg::distribute(void)
{
  int n = srcPatchID.size();
  int m = 0;
  for ( int i = 0; i < n; ++i )
  {
    MigrateAtomsMsg *msg = new MigrateAtomsMsg;
    msg->fromNodeID = fromNodeID;
    msg->srcPatchID = srcPatchID[i];
    msg->destPatchID = destPatchID[i];
    int l = numAtoms[i];
    if ( l )
    {
      DebugM(3,"Distributing " << l << " atoms to patch " << msg->destPatchID << "\n");
      msg->migrationList = new MigrationList(l);
      for ( int j = 0; j < l; ++j ) (*msg->migrationList)[j] = migrationList[m+j];
      m += l;
    }
    else
    {
      msg->migrationList = 0;
    }
    PatchMap::Object()->homePatch(msg->destPatchID)->depositMigration(msg);
  }
}

void * MigrateAtomsCombinedMsg::pack (MigrateAtomsCombinedMsg *m) {
  int n = m->srcPatchID.size();
  int l = sizeof(NodeID)			// fromNodeID
	+ sizeof(int)				// totalAtoms
	+ sizeof(int)				// n
	+ sizeof(PatchID) * n			// srcPatchID
	+ sizeof(PatchID) * n			// destPatchID
	+ sizeof(int) * n			// numAtoms
	+ sizeof(MigrationElem) * m->totalAtoms	// migrationList
  ;

  char *buffer;
  char *b = buffer = (char*)CkAllocBuffer(m,l);

  *((NodeID*)b) = m->fromNodeID;	b += sizeof(NodeID);
  *((int*)b) = m->totalAtoms;	b += sizeof(int);
  *((int*)b) = n;		b += sizeof(int);

  memcpy(b,(void*)&(m->srcPatchID[0]), sizeof(PatchID) * n);
  b += sizeof(PatchID) * n;

  memcpy(b,(void*)&(m->destPatchID[0]), sizeof(PatchID) * n);
  b += sizeof(PatchID) * n;

  memcpy(b,(void*)&(m->numAtoms[0]), sizeof(int) * n);
  b += sizeof(int) * n;

  memcpy(b,(void*)&(m->migrationList[0]),sizeof(MigrationElem) * m->totalAtoms);
  b += sizeof(MigrationElem) * m->totalAtoms;

  delete m;
  return buffer;
}

MigrateAtomsCombinedMsg* MigrateAtomsCombinedMsg::unpack(void *ptr) {
  void *_ptr = CkAllocBuffer(ptr, sizeof(MigrateAtomsCombinedMsg));
  MigrateAtomsCombinedMsg* m = new (_ptr) MigrateAtomsCombinedMsg;
  char *b = (char*)ptr;

  m->fromNodeID = *((NodeID*)b);	b += sizeof(NodeID);
  m->totalAtoms = *((int*)b);	b += sizeof(int);
  int n = *((int*)b);		b += sizeof(int);

  DebugM(3,"Unpacking MigrateAtomsCombinedMsg with " << n << " messages.\n");

  m->srcPatchID.resize(n);
  memcpy((void*)&(m->srcPatchID[0]), b, sizeof(PatchID) * n);
  b += sizeof(PatchID) * n;

  m->destPatchID.resize(n);
  memcpy((void*)&(m->destPatchID[0]), b, sizeof(PatchID) * n);
  b += sizeof(PatchID) * n;

  m->numAtoms.resize(n);
  memcpy((void*)&(m->numAtoms[0]), b, sizeof(int) * n);
  b += sizeof(int) * n;

  m->migrationList.resize(m->totalAtoms);
  memcpy((void*)&(m->migrationList[0]), b, sizeof(MigrationElem)*m->totalAtoms);
  b += sizeof(MigrationElem) * m->totalAtoms;
  CkFreeMsg(ptr);
  return m;
}



/***************************************************************************
 * RCS INFORMATION:
 *
 *	$RCSfile: MigrateAtomsMsg.C,v $
 *	$Author: brunner $	$Locker:  $		$State: Exp $
 *	$Revision: 1.9 $	$Date: 1999/05/11 23:56:35 $
 *
 * REVISION HISTORY:
 *
 * $Log: MigrateAtomsMsg.C,v $
 * Revision 1.9  1999/05/11 23:56:35  brunner
 * Changes for new charm version
 *
 * Revision 1.8  1998/10/24 19:57:39  jim
 * Eliminated warnings generated by g++ -Wall.
 *
 * Revision 1.7  1998/03/03 23:05:16  brunner
 * Changed include files for new simplified Charm++ include file structure.
 *
 * Revision 1.6  1997/04/11 16:54:28  jim
 * Fixed bug calling size() on possibly null pointer.
 *
 * Revision 1.5  1997/04/11 06:03:22  jim
 * Message combining implemented for atom migration.
 *
 * Revision 1.4  1997/04/10 22:29:12  jim
 * First steps towards combining atom migration messages.
 *
 * Revision 1.3  1997/03/06 22:06:04  ari
 * Removed Compute.ci
 * Comments added - more code cleaning
 *
 * Revision 1.2  1997/02/26 16:53:11  ari
 * Cleaning and debuging for memory leaks.
 * Adding comments.
 * Removed some dead code due to use of Quiescense detection.
 *
 * Revision 1.1  1997/02/17 23:47:01  ari
 * Added files for cleaning up atom migration code
 *
 ***************************************************************************/
