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

#ifndef MAIN_H
#define MAIN_H

#include "charm++.h"

#include "BOCgroup.h"
#include "NamdTypes.h"

#include "main.decl.h"

// Message to send our per processor BOC's list of groupIDs of 
// all other BOC's
class GroupInitMsg : public CMessage_GroupInitMsg
{
public:
  BOCgroup group;
};

class SlaveInitMsg : public GroupInitMsg
{
public:
  CkChareID master;
};

class Compute;

// For Compute objects to enqueue themselves when ready to compute
class LocalWorkMsg : public CMessage_LocalWorkMsg
{
public:
  Compute *compute;
};

#endif /* MAIN_H */

