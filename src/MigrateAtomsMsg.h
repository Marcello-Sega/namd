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
