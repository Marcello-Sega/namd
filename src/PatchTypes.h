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
  int doFullElectrostatics;

private:
  int spacer;  // Use this to keep byte-aligned for now.  -JCP
};

class Results
{
public:
  enum { normal, slow, maxNumForces };
  Force *f[maxNumForces];
};

#endif
