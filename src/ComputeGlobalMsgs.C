/***************************************************************************/
/*         (C) Copyright 1996,1997 The Board of Trustees of the            */
/*                          University of Illinois                         */
/*                           All Rights Reserved                           */
/***************************************************************************/
/***************************************************************************
 * DESCRIPTION:  Messages needed for ComputeGlobal operation.
 *		
 ***************************************************************************/

static char ident[] = "@(#)$Header: /home/cvs/namd/cvsroot/namd2/src/ComputeGlobalMsgs.C,v 1.2 1998/02/16 00:24:38 jim Exp $";

#include "ComputeGlobalMsgs.h"

//#define DEBUGM
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
  int gsize = gdef.size();
  *length = 2 * sizeof(int) + (size+gsize) * sizeof(AtomID);

  char *buffer;
  char *b = buffer = (char*)new_packbuffer(this,*length);

  memcpy(b, &size, sizeof(int)); b += sizeof(int);
  memcpy(b, aid.begin(), size*sizeof(AtomID)); b += size*sizeof(AtomID);
  memcpy(b, &gsize, sizeof(int)); b += sizeof(int);
  memcpy(b, gdef.begin(), gsize*sizeof(AtomID)); b += gsize*sizeof(AtomID);

  this->~ComputeGlobalConfigMsg();
  return buffer;
}

void ComputeGlobalConfigMsg::unpack (void *in) {
  new((void*)this) ComputeGlobalConfigMsg;
  char *b = (char*)in;

  int size, gsize;
  memcpy(&size, b, sizeof(int)); b += sizeof(int);
  aid.resize(size);
  memcpy(aid.begin(), b, size*sizeof(AtomID)); b += size*sizeof(AtomID);
  memcpy(&gsize, b, sizeof(int)); b += sizeof(int);
  gdef.resize(gsize);
  memcpy(gdef.begin(), b, gsize*sizeof(AtomID)); b += gsize*sizeof(AtomID);

  // DO NOT delete void *in - this is done by Charm
}


// DATA MESSAGE

ComputeGlobalDataMsg::ComputeGlobalDataMsg(void) { 
}

ComputeGlobalDataMsg::~ComputeGlobalDataMsg(void) { 
}

void * ComputeGlobalDataMsg::pack (int *length) {
  int size = aid.size();
  int gsize = gcom.size();
  *length = 2 * sizeof(int) + size * (sizeof(AtomID)+sizeof(Position))
			    + gsize * sizeof(Position);

  char *buffer;
  char *b = buffer = (char*)new_packbuffer(this,*length);

  memcpy(b, &size, sizeof(int)); b += sizeof(int);
  memcpy(b, aid.begin(), size*sizeof(AtomID)); b += size*sizeof(AtomID);
  memcpy(b, p.begin(), size*sizeof(Position)); b += size*sizeof(Position);
  memcpy(b, &gsize, sizeof(int)); b += sizeof(int);
  memcpy(b, gcom.begin(), gsize*sizeof(Position)); b += gsize*sizeof(Position);

  this->~ComputeGlobalDataMsg();
  return buffer;
}

void ComputeGlobalDataMsg::unpack (void *in) {
  new((void*)this) ComputeGlobalDataMsg;
  char *b = (char*)in;

  int size, gsize;
  memcpy(&size, b, sizeof(int)); b += sizeof(int);
  aid.resize(size);
  memcpy(aid.begin(), b, size*sizeof(AtomID)); b += size*sizeof(AtomID);
  p.resize(size);
  memcpy(p.begin(), b, size*sizeof(Position)); b += size*sizeof(Position);
  memcpy(&gsize, b, sizeof(int)); b += sizeof(int);
  gcom.resize(gsize);
  memcpy(gcom.begin(), b, gsize*sizeof(Position)); b += gsize*sizeof(Position);

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
  int gsize = gforce.size();
  *length = 3*sizeof(int) + size * (sizeof(AtomID)+sizeof(Force))
			  + gsize * sizeof(Force);
  if ( reconfig ) *length += sizeof(int) + newaid.size() * sizeof(AtomID) +
		             sizeof(int) + newgdef.size() * sizeof(AtomID);

  char *buffer;
  char *b = buffer = (char*)new_packbuffer(this,*length);

  memcpy(b, &size, sizeof(int)); b += sizeof(int);
  memcpy(b, aid.begin(), size*sizeof(AtomID)); b += size*sizeof(AtomID);
  memcpy(b, f.begin(), size*sizeof(Force)); b += size*sizeof(Force);
  memcpy(b, &gsize, sizeof(int)); b += sizeof(int);
  memcpy(b, gforce.begin(), gsize*sizeof(Force)); b += gsize*sizeof(Force);
  memcpy(b, &reconfig, sizeof(int)); b += sizeof(int);
  if ( reconfig ) {
    size = newaid.size();
    memcpy(b, &size, sizeof(int)); b += sizeof(int);
    memcpy(b, newaid.begin(), size*sizeof(AtomID)); b += size*sizeof(AtomID);
    gsize = newgdef.size();
    memcpy(b, &gsize, sizeof(int)); b += sizeof(int);
    memcpy(b, newgdef.begin(), gsize*sizeof(AtomID)); b += gsize*sizeof(AtomID);
  }

  this->~ComputeGlobalResultsMsg();
  return buffer;
}

void ComputeGlobalResultsMsg::unpack (void *in) {
  new((void*)this) ComputeGlobalResultsMsg;
  char *b = (char*)in;

  int size, gsize;
  memcpy(&size, b, sizeof(int)); b += sizeof(int);
  aid.resize(size);
  memcpy(aid.begin(), b, size*sizeof(AtomID)); b += size*sizeof(AtomID);
  f.resize(size);
  memcpy(f.begin(), b, size*sizeof(Force)); b += size*sizeof(Force);
  memcpy(&gsize, b, sizeof(int)); b += sizeof(int);
  gforce.resize(gsize);
  memcpy(gforce.begin(), b, gsize*sizeof(Force)); b += gsize*sizeof(Force);
  memcpy(&reconfig, b, sizeof(int)); b += sizeof(int);
  if ( reconfig ) {
    memcpy(&size, b, sizeof(int)); b += sizeof(int);
    newaid.resize(size);
    memcpy(newaid.begin(), b, size*sizeof(AtomID)); b += size*sizeof(AtomID);
    memcpy(&gsize, b, sizeof(int)); b += sizeof(int);
    newgdef.resize(gsize);
    memcpy(newgdef.begin(), b, gsize*sizeof(AtomID)); b += gsize*sizeof(AtomID);
  }

  // DO NOT delete void *in - this is done by Charm
}



/***************************************************************************
 * RCS INFORMATION:
 *
 *	$RCSfile: ComputeGlobalMsgs.C,v $
 *	$Author: jim $	$Locker:  $		$State: Exp $
 *	$Revision: 1.2 $	$Date: 1998/02/16 00:24:38 $
 *
 * REVISION HISTORY:
 *
 * $Log: ComputeGlobalMsgs.C,v $
 * Revision 1.2  1998/02/16 00:24:38  jim
 * Added atom group centers of mass to Tcl interface.
 *
 * Revision 1.1  1997/12/19 23:48:47  jim
 * Added Tcl interface for calculating forces.
 *
 *
 ***************************************************************************/
