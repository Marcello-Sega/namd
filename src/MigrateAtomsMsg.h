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

#ifndef MIGRATEATOMSMSG_H
#define MIGRATEATOMSMSG_H

#include "ckdefs.h"
#include "chare.h"
#include "c++interface.h"

#include "NamdTypes.h"
#include "Templates/SortedArray.h"
#include "Migration.h"

// Message which stores list of atoms and their data
// which are to be migrated from one patch to another.
// This message does not contain information that will change asynchronously
// It does not need to be prepacked
class MigrateAtomsMsg : public comm_object {
public:
  NodeID  fromNodeID;
  PatchID srcPatchID;
  PatchID destPatchID;
  MigrationList *migrationList;

  MigrateAtomsMsg(void);
  ~MigrateAtomsMsg(void);

  MigrateAtomsMsg(PatchID source, PatchID destination, MigrationList *m);

  void * operator new(size_t s, int i) {return comm_object::operator new(s,i);}
  void * operator new(size_t s) { return comm_object::operator new(s); }
  void * operator new(size_t, void *ptr) { return ptr; }

  // pack and unpack functions
  void * pack (int *length);
  void unpack (void *in);
};

#endif
