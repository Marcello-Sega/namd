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

#include "charm++.h"

#include "NamdTypes.h"
#include "SortedArray.h"
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

class MigrateAtomsCombinedMsg : public comm_object
{
public:
  NodeID fromNodeID;
  ResizeArray<PatchID> srcPatchID;
  ResizeArray<PatchID> destPatchID;
  ResizeArray<int> numAtoms;
  int totalAtoms;
  MigrationList migrationList;

  MigrateAtomsCombinedMsg(void);
  ~MigrateAtomsCombinedMsg(void) { };

  void add(PatchID source, PatchID destination, MigrationList *m);
  void distribute(void);

  // Standard new overload for comm_object new
  void * operator new(size_t s, int i) {return comm_object::operator new(s,i);}
  void * operator new(size_t s) { return comm_object::operator new(s); }
  void * operator new(size_t, void *ptr) { return ptr; }

  // pack and unpack functions
  void * pack (int *length);
  void unpack (void *in);
};

#endif


/***************************************************************************
 * RCS INFORMATION:
 *
 *	$RCSfile $
 *	$Author $	$Locker:  $		$State: Exp $
 *	$Revision: 1.8 $	$Date: 1998/03/03 23:05:16 $
 *
 ***************************************************************************
 * REVISION HISTORY:
 *
 * $Log: MigrateAtomsMsg.h,v $
 * Revision 1.8  1998/03/03 23:05:16  brunner
 * Changed include files for new simplified Charm++ include file structure.
 *
 * Revision 1.7  1997/12/26 23:10:50  milind
 * Made namd2 to compile, link and run under linux. Merged Templates and src
 * directoriies, and removed separate definition and declaration files for
 * templates.
 *
 * Revision 1.6  1997/04/11 06:03:24  jim
 * Message combining implemented for atom migration.
 *
 * Revision 1.5  1997/04/10 22:29:14  jim
 * First steps towards combining atom migration messages.
 *
 * Revision 1.4  1997/03/19 11:54:30  ari
 * Add Broadcast mechanism.
 * Fixed RCS Log entries on files that did not have Log entries.
 * Added some register variables to Molecule and ComputeNonbondedExcl.C
 *
 ***************************************************************************/
