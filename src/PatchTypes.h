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
#include "Lattice.h"

class Flags
{
public:
  int step;			// timestep number reported to user
  				// Same number may appear multiple times!
  int doNonbonded;
  int doFullElectrostatics;
  int doMolly;
  int submitLoadStats;
  int maxForceUsed;		// may ignore slower force classes
  int maxForceMerged;		// add this and faster to normal

  Lattice lattice;		// rather than shipping around separately

};

class Results
{
public:
  enum { normal=0, nbond=1, slow=2, maxNumForces=3 };
  Force *f[maxNumForces];
};

#endif

