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
  int doFullElectrostatics;

private:
  // int spacer;  // Use this to keep byte-aligned for now.  -JCP
};

class Results
{
public:
  enum { normal, slow, maxNumForces };
  Force *f[maxNumForces];
};

#endif


/***************************************************************************
 * RCS INFORMATION:
 *
 *	$RCSfile $
 *	$Author $	$Locker:  $		$State: Exp $
 *	$Revision: 1.4 $	$Date: 1997/03/19 11:54:48 $
 *
 ***************************************************************************
 * REVISION HISTORY:
 *
 * $Log: PatchTypes.h,v $
 * Revision 1.4  1997/03/19 11:54:48  ari
 * Add Broadcast mechanism.
 * Fixed RCS Log entries on files that did not have Log entries.
 * Added some register variables to Molecule and ComputeNonbondedExcl.C
 *
 ***************************************************************************/
