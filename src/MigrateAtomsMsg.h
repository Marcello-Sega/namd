//-*-c++-*-
/***************************************************************************/
/*         (C) Copyright 1996,1997 The Board of Trustees of the            */
/*                          University of Illinois                         */
/*                           All Rights Reserved                           */
/***************************************************************************/
/***************************************************************************
 * DESCRIPTION: Migration messages are sent from HomePatch to HomePatch
 *              with lists of atoms and atom information (if any) that
 *              need to be migrated.  A message must be sent from a
 *              neighbor even if null so that the HomePatch knows
 *              what atoms it will have before commencing a positionsReady()
 *              to its Computes.
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

  // Standard new overload for comm_object new
  void * operator new(size_t s, int i) {return comm_object::operator new(s,i);}
  void * operator new(size_t s) { return comm_object::operator new(s); }
  void * operator new(size_t, void *ptr) { return ptr; }

  // pack and unpack functions
  void * pack (int *length);
  void unpack (void *in);
};

#endif
