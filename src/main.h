/**
***  Copyright (c) 1995, 1996, 1997, 1998, 1999, 2000 by
***  The Board of Trustees of the University of Illinois.
***  All rights reserved.
**/

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

