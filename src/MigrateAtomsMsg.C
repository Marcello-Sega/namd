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

static char ident[] = "@(#)$Header: /home/cvs/namd/cvsroot/namd2/src/MigrateAtomsMsg.C,v 1.7 1998/03/03 23:05:16 brunner Exp $";

#include "charm++.h"

#include "Migration.h"
#include "MigrateAtomsMsg.h"
#include "NamdTypes.h"
// #define DEBUGM
#define MIN_DEBUG_LEVEL 3
#include "Debug.h"

#include "PatchMgr.top.h"
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
    fromNodeID = CMyPe();
}

void * MigrateAtomsMsg::pack (int *length) {
  *length = sizeof(NodeID) + sizeof(PatchID) + sizeof(PatchID) + sizeof(int) 
    + (migrationList != NULL 
	   ? migrationList->size() * sizeof(MigrationElem) 
	   : 0 );

  // This is what I want
  char *buffer;
  char *b = buffer = (char*)new_packbuffer(this,*length);

  *((NodeID*)b) = fromNodeID; b += sizeof(NodeID);
  *((PatchID*)b) = srcPatchID; b += sizeof(PatchID);
  *((PatchID*)b) = destPatchID; b += sizeof(PatchID);

  if (migrationList != NULL) {
    *((int*)b) = migrationList->size(); b += sizeof(int);
    for ( int i = 0; i < migrationList->size(); i++ ) {
      *((MigrationElem*)b) = (*migrationList)[i]; b += sizeof(MigrationElem);
    }
  }
  else {
    *((int*)b) = 0; b += sizeof(int);
  }

  // Delete of new'd item from HomePatch::doAtomMigration()
  delete migrationList;
  migrationList = 0;

  this->~MigrateAtomsMsg();
  return buffer;
}

void MigrateAtomsMsg::unpack (void *in) {
  new((void*)this) MigrateAtomsMsg;
  char *b = (char*)in;
  fromNodeID = *((NodeID*)b); b += sizeof(NodeID);
  srcPatchID = *((PatchID*)b); b += sizeof(PatchID);
  destPatchID = *((PatchID*)b); b += sizeof(PatchID);
  int size = *((int*)b); b += sizeof(int);
  if (size != 0) {
    migrationList = new MigrationList();
    migrationList->resize(size);
    for ( int i = 0; i < size; i++ ) {
      (*migrationList)[i] = *((MigrationElem*)b); b += sizeof(MigrationElem);
    }
  }
  else {
    migrationList = NULL;
  }
  // DO NOT delete void *in - this is done by Charm
}

MigrateAtomsCombinedMsg::MigrateAtomsCombinedMsg(void)
{
  fromNodeID = CMyPe();
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
    MigrateAtomsMsg *msg = new (MsgIndex(MigrateAtomsMsg)) MigrateAtomsMsg;
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

void * MigrateAtomsCombinedMsg::pack (int *length) {
  int n = srcPatchID.size();
  int l = sizeof(NodeID)			// fromNodeID
	+ sizeof(int)				// totalAtoms
	+ sizeof(int)				// n
	+ sizeof(PatchID) * n			// srcPatchID
	+ sizeof(PatchID) * n			// destPatchID
	+ sizeof(int) * n			// numAtoms
	+ sizeof(MigrationElem) * totalAtoms	// migrationList
  ; *length = l;

  char *buffer;
  char *b = buffer = (char*)new_packbuffer(this,*length);

  *((NodeID*)b) = fromNodeID;	b += sizeof(NodeID);
  *((int*)b) = totalAtoms;	b += sizeof(int);
  *((int*)b) = n;		b += sizeof(int);

  memcpy(b,(void*)&srcPatchID[0], sizeof(PatchID) * n);
  b += sizeof(PatchID) * n;

  memcpy(b,(void*)&destPatchID[0], sizeof(PatchID) * n);
  b += sizeof(PatchID) * n;

  memcpy(b,(void*)&numAtoms[0], sizeof(int) * n);
  b += sizeof(int) * n;

  memcpy(b,(void*)&migrationList[0], sizeof(MigrationElem) * totalAtoms);
  b += sizeof(MigrationElem) * totalAtoms;

  this->~MigrateAtomsCombinedMsg();
  return buffer;
}

void MigrateAtomsCombinedMsg::unpack (void *in) {
  new((void*)this) MigrateAtomsCombinedMsg;
  char *b = (char*)in;

  fromNodeID = *((NodeID*)b);	b += sizeof(NodeID);
  totalAtoms = *((int*)b);	b += sizeof(int);
  int n = *((int*)b);		b += sizeof(int);

  DebugM(3,"Unpacking MigrateAtomsCombinedMsg with " << n << " messages.\n");

  srcPatchID.resize(n);
  memcpy((void*)&srcPatchID[0], b, sizeof(PatchID) * n);
  b += sizeof(PatchID) * n;

  destPatchID.resize(n);
  memcpy((void*)&destPatchID[0], b, sizeof(PatchID) * n);
  b += sizeof(PatchID) * n;

  numAtoms.resize(n);
  memcpy((void*)&numAtoms[0], b, sizeof(int) * n);
  b += sizeof(int) * n;

  migrationList.resize(totalAtoms);
  memcpy((void*)&migrationList[0], b, sizeof(MigrationElem) * totalAtoms);
  b += sizeof(MigrationElem) * totalAtoms;

  // DO NOT delete void *in - this is done by Charm
}



/***************************************************************************
 * RCS INFORMATION:
 *
 *	$RCSfile: MigrateAtomsMsg.C,v $
 *	$Author: brunner $	$Locker:  $		$State: Exp $
 *	$Revision: 1.7 $	$Date: 1998/03/03 23:05:16 $
 *
 * REVISION HISTORY:
 *
 * $Log: MigrateAtomsMsg.C,v $
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
