//-*-c++-*-
/***************************************************************************/
/*                                                                         */
/*              (C) Copyright 1996 The Board of Trustees of the            */
/*                          University of Illinois                         */
/*                           All Rights Reserved                           */
/*                                                                         */
/***************************************************************************/

/***************************************************************************
 * DESCRIPTION:
 *
 ***************************************************************************/

#ifndef PATCHTYPES_H
#define PATCHTYPES_H

#include "NamdTypes.h"

class Flags
{
public:
  int seq;			// sequence number
  int doNonbonded;
  int doFullElectrostatics;
  int submitLoadStats;

private:
//int spacer;  // Use this to keep byte-aligned for now.  -JCP
               // Actually double-word aligned, I think -RKB
};

class Results
{
public:
  enum { normal, nbond, slow, maxNumForces };
  Force *f[maxNumForces];
};

#endif


/***************************************************************************
 * RCS INFORMATION:
 *
 *	$RCSfile $
 *	$Author $	$Locker:  $		$State: Exp $
 *	$Revision: 1.6 $	$Date: 1997/03/27 20:25:50 $
 *
 ***************************************************************************
 * REVISION HISTORY:
 *
 * $Log: PatchTypes.h,v $
 * Revision 1.6  1997/03/27 20:25:50  brunner
 * Changes for LdbCoordinator, the load balance control BOC
 *
 * Revision 1.5  1997/03/25 23:01:00  jim
 * Added nonbondedFrequency parameter and multiple time-stepping
 *
 * Revision 1.4  1997/03/19 11:54:48  ari
 * Add Broadcast mechanism.
 * Fixed RCS Log entries on files that did not have Log entries.
 * Added some register variables to Molecule and ComputeNonbondedExcl.C
 *
 ***************************************************************************/
