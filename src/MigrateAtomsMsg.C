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

static char ident[] = "@(#)$Header: /home/cvs/namd/cvsroot/namd2/src/MigrateAtomsMsg.C,v 1.3 1997/03/06 22:06:04 ari Exp $";

#include "ckdefs.h"
#include "chare.h"
#include "c++interface.h"
#include "Migration.h"
#include "MigrateAtomsMsg.h"
#include "NamdTypes.h"
// #define DEBUGM
#define MIN_DEBUG_LEVEL 3
#include "Debug.h"

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

/***************************************************************************
 * RCS INFORMATION:
 *
 *	$RCSfile: MigrateAtomsMsg.C,v $
 *	$Author: ari $	$Locker:  $		$State: Exp $
 *	$Revision: 1.3 $	$Date: 1997/03/06 22:06:04 $
 *
 * REVISION HISTORY:
 *
 * $Log: MigrateAtomsMsg.C,v $
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
