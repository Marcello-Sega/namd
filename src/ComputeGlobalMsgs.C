/***************************************************************************/
/*         (C) Copyright 1996,1997 The Board of Trustees of the            */
/*                          University of Illinois                         */
/*                           All Rights Reserved                           */
/***************************************************************************/
/***************************************************************************
 * DESCRIPTION:  Messages needed for ComputeGlobal operation.
 *		
 ***************************************************************************/

static char ident[] = "@(#)$Header: /home/cvs/namd/cvsroot/namd2/src/ComputeGlobalMsgs.C,v 1.1 1997/12/19 23:48:47 jim Exp $";

#include "ComputeGlobalMsgs.h"

// #define DEBUGM
#define MIN_DEBUG_LEVEL 3
#include "Debug.h"

#include "ComputeMgr.top.h"


// CONFIG MESSAGE

ComputeGlobalConfigMsg::ComputeGlobalConfigMsg(void) { 
}

ComputeGlobalConfigMsg::~ComputeGlobalConfigMsg(void) { 
}

void * ComputeGlobalConfigMsg::pack (int *length) {
  int size = aid.size();
  *length = sizeof(int) + size * sizeof(AtomID);

  char *buffer;
  char *b = buffer = (char*)new_packbuffer(this,*length);

  memcpy(b, &size, sizeof(int)); b += sizeof(int);
  memcpy(b, aid.begin(), size*sizeof(AtomID)); b += size*sizeof(AtomID);

  this->~ComputeGlobalConfigMsg();
  return buffer;
}

void ComputeGlobalConfigMsg::unpack (void *in) {
  new((void*)this) ComputeGlobalConfigMsg;
  char *b = (char*)in;

  int size;
  memcpy(&size, b, sizeof(int)); b += sizeof(int);
  aid.resize(size);
  memcpy(aid.begin(), b, size*sizeof(AtomID)); b += size*sizeof(AtomID);

  // DO NOT delete void *in - this is done by Charm
}


// DATA MESSAGE

ComputeGlobalDataMsg::ComputeGlobalDataMsg(void) { 
}

ComputeGlobalDataMsg::~ComputeGlobalDataMsg(void) { 
}

void * ComputeGlobalDataMsg::pack (int *length) {
  int size = aid.size();
  *length = sizeof(int) + size * (sizeof(AtomID)+sizeof(Position));

  char *buffer;
  char *b = buffer = (char*)new_packbuffer(this,*length);

  memcpy(b, &size, sizeof(int)); b += sizeof(int);
  memcpy(b, aid.begin(), size*sizeof(AtomID)); b += size*sizeof(AtomID);
  memcpy(b, p.begin(), size*sizeof(Position)); b += size*sizeof(Position);

  this->~ComputeGlobalDataMsg();
  return buffer;
}

void ComputeGlobalDataMsg::unpack (void *in) {
  new((void*)this) ComputeGlobalDataMsg;
  char *b = (char*)in;

  int size;
  memcpy(&size, b, sizeof(int)); b += sizeof(int);
  aid.resize(size);
  memcpy(aid.begin(), b, size*sizeof(AtomID)); b += size*sizeof(AtomID);
  p.resize(size);
  memcpy(p.begin(), b, size*sizeof(Position)); b += size*sizeof(Position);

  // DO NOT delete void *in - this is done by Charm
}


// RESULTS MESSAGE

ComputeGlobalResultsMsg::ComputeGlobalResultsMsg(void) { 
  reconfig = 0;
}

ComputeGlobalResultsMsg::~ComputeGlobalResultsMsg(void) { 
}

void * ComputeGlobalResultsMsg::pack (int *length) {
  int size = aid.size();
  *length = 2*sizeof(int) + size * (sizeof(AtomID)+sizeof(Force));
  if ( reconfig ) *length += sizeof(int) + newaid.size() * sizeof(AtomID);

  char *buffer;
  char *b = buffer = (char*)new_packbuffer(this,*length);

  memcpy(b, &size, sizeof(int)); b += sizeof(int);
  memcpy(b, aid.begin(), size*sizeof(AtomID)); b += size*sizeof(AtomID);
  memcpy(b, f.begin(), size*sizeof(Force)); b += size*sizeof(Force);
  memcpy(b, &reconfig, sizeof(int)); b += sizeof(int);
  if ( reconfig ) {
    size = newaid.size();
    memcpy(b, &size, sizeof(int)); b += sizeof(int);
    memcpy(b, newaid.begin(), size*sizeof(AtomID)); b += size*sizeof(AtomID);
  }

  this->~ComputeGlobalResultsMsg();
  return buffer;
}

void ComputeGlobalResultsMsg::unpack (void *in) {
  new((void*)this) ComputeGlobalResultsMsg;
  char *b = (char*)in;

  int size;
  memcpy(&size, b, sizeof(int)); b += sizeof(int);
  aid.resize(size);
  memcpy(aid.begin(), b, size*sizeof(AtomID)); b += size*sizeof(AtomID);
  f.resize(size);
  memcpy(f.begin(), b, size*sizeof(Force)); b += size*sizeof(Force);
  memcpy(&reconfig, b, sizeof(int)); b += sizeof(int);
  if ( reconfig ) {
    memcpy(&size, b, sizeof(int)); b += sizeof(int);
    newaid.resize(size);
    memcpy(newaid.begin(), b, size*sizeof(AtomID)); b += size*sizeof(AtomID);
  }

  // DO NOT delete void *in - this is done by Charm
}



/***************************************************************************
 * RCS INFORMATION:
 *
 *	$RCSfile: ComputeGlobalMsgs.C,v $
 *	$Author: jim $	$Locker:  $		$State: Exp $
 *	$Revision: 1.1 $	$Date: 1997/12/19 23:48:47 $
 *
 * REVISION HISTORY:
 *
 * $Log: ComputeGlobalMsgs.C,v $
 * Revision 1.1  1997/12/19 23:48:47  jim
 * Added Tcl interface for calculating forces.
 *
 *
 ***************************************************************************/
