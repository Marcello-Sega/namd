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

class Flags
{
public:
  int doFullElectrostatics;

private:
  int spacer;  // Use this to keep byte-aligned for now.  -JCP
};

#endif
