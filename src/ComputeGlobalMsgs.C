/***************************************************************************/
/*         (C) Copyright 1996,1997 The Board of Trustees of the            */
/*                          University of Illinois                         */
/*                           All Rights Reserved                           */
/***************************************************************************/
/***************************************************************************
 * DESCRIPTION:  Messages needed for ComputeGlobal operation.
 *		
 ***************************************************************************/

#include "ComputeGlobalMsgs.h"

//#define DEBUGM
#define MIN_DEBUG_LEVEL 3
#include "Debug.h"

#include "ComputeMgr.decl.h"


// CONFIG MESSAGE

ComputeGlobalConfigMsg::ComputeGlobalConfigMsg(void) { 
}

ComputeGlobalConfigMsg::~ComputeGlobalConfigMsg(void) { 
}

void* ComputeGlobalConfigMsg::pack (ComputeGlobalConfigMsg *msg) {
  int size = msg->aid.size();
  int gsize = msg->gdef.size();
  int length = 2 * sizeof(int) + (size+gsize) * sizeof(AtomID);

  char *buffer;
  char *b = buffer = (char*)CkAllocBuffer(msg,length);

  memcpy(b, &size, sizeof(int)); b += sizeof(int);
  memcpy(b, msg->aid.begin(), size*sizeof(AtomID)); b += size*sizeof(AtomID);
  memcpy(b, &gsize, sizeof(int)); b += sizeof(int);
  memcpy(b, msg->gdef.begin(), gsize*sizeof(AtomID)); b += gsize*sizeof(AtomID);

  delete msg;
  return buffer;
}

ComputeGlobalConfigMsg* ComputeGlobalConfigMsg::unpack (void *ptr) {
  void *_ptr = CkAllocBuffer(ptr, sizeof(ComputeGlobalConfigMsg));
  ComputeGlobalConfigMsg* m = new (_ptr) ComputeGlobalConfigMsg;
  char *b = (char*)ptr;

  int size, gsize;
  memcpy(&size, b, sizeof(int)); b += sizeof(int);
  m->aid.resize(size);
  memcpy(m->aid.begin(), b, size*sizeof(AtomID)); b += size*sizeof(AtomID);
  memcpy(&gsize, b, sizeof(int)); b += sizeof(int);
  m->gdef.resize(gsize);
  memcpy(m->gdef.begin(), b, gsize*sizeof(AtomID)); b += gsize*sizeof(AtomID);
  CkFreeMsg(ptr);
  return m;
}


// DATA MESSAGE

ComputeGlobalDataMsg::ComputeGlobalDataMsg(void) { 
}

ComputeGlobalDataMsg::~ComputeGlobalDataMsg(void) { 
}

void * ComputeGlobalDataMsg::pack (ComputeGlobalDataMsg *m) {
  int size = m->aid.size();
  int gsize = m->gcom.size();
  int length = 2 * sizeof(int) + size * (sizeof(AtomID)+sizeof(Position))
			    + gsize * sizeof(Position);

  char *buffer;
  char *b = buffer = (char*)CkAllocBuffer(m,length);

  memcpy(b, &size, sizeof(int)); b += sizeof(int);
  memcpy(b, m->aid.begin(), size*sizeof(AtomID)); b += size*sizeof(AtomID);
  memcpy(b, m->p.begin(), size*sizeof(Position)); b += size*sizeof(Position);
  memcpy(b, &gsize, sizeof(int)); b += sizeof(int);
  memcpy(b, m->gcom.begin(), gsize*sizeof(Position)); b += gsize*sizeof(Position);

  delete m;
  return buffer;
}

ComputeGlobalDataMsg* ComputeGlobalDataMsg::unpack (void *ptr) {
  void *_ptr = CkAllocBuffer(ptr, sizeof(ComputeGlobalDataMsg));
  ComputeGlobalDataMsg* m = new (_ptr) ComputeGlobalDataMsg;
  char *b = (char*)ptr;

  int size, gsize;
  memcpy(&size, b, sizeof(int)); b += sizeof(int);
  m->aid.resize(size);
  memcpy(m->aid.begin(), b, size*sizeof(AtomID)); b += size*sizeof(AtomID);
  m->p.resize(size);
  memcpy(m->p.begin(), b, size*sizeof(Position)); b += size*sizeof(Position);
  memcpy(&gsize, b, sizeof(int)); b += sizeof(int);
  m->gcom.resize(gsize);
  memcpy(m->gcom.begin(), b, gsize*sizeof(Position)); b += gsize*sizeof(Position);
  CkFreeMsg(ptr);
  return m;
}


// RESULTS MESSAGE

ComputeGlobalResultsMsg::ComputeGlobalResultsMsg(void) { 
  reconfig = 0;
}

ComputeGlobalResultsMsg::~ComputeGlobalResultsMsg(void) { 
}

void* ComputeGlobalResultsMsg::pack (ComputeGlobalResultsMsg *msg) {
  int size = msg->aid.size();
  int gsize = msg->gforce.size();
  int length = 3*sizeof(int) + size * (sizeof(AtomID)+sizeof(Force))
			  + gsize * sizeof(Force);
  if ( msg->reconfig ) length += sizeof(int) + 
                             msg->newaid.size() * sizeof(AtomID) +
		             sizeof(int) + msg->newgdef.size() * sizeof(AtomID);

  char *buffer;
  char *b = buffer = (char*)CkAllocBuffer(msg,length);

  memcpy(b, &size, sizeof(int)); b += sizeof(int);
  memcpy(b, msg->aid.begin(), size*sizeof(AtomID)); b += size*sizeof(AtomID);
  memcpy(b, msg->f.begin(), size*sizeof(Force)); b += size*sizeof(Force);
  memcpy(b, &gsize, sizeof(int)); b += sizeof(int);
  memcpy(b, msg->gforce.begin(), gsize*sizeof(Force)); b += gsize*sizeof(Force);
  memcpy(b, &(msg->reconfig), sizeof(int)); b += sizeof(int);
  if ( msg->reconfig ) {
    size = msg->newaid.size();
    memcpy(b, &size, sizeof(int)); b += sizeof(int);
    memcpy(b, msg->newaid.begin(), size*sizeof(AtomID)); b += size*sizeof(AtomID);
    gsize = msg->newgdef.size();
    memcpy(b, &gsize, sizeof(int)); b += sizeof(int);
    memcpy(b, msg->newgdef.begin(), gsize*sizeof(AtomID)); b += gsize*sizeof(AtomID);
  }
  delete msg;
  return buffer;
}

ComputeGlobalResultsMsg* ComputeGlobalResultsMsg::unpack(void *ptr) {
  void *_ptr = CkAllocBuffer(ptr, sizeof(ComputeGlobalResultsMsg));
  ComputeGlobalResultsMsg* m = new (_ptr) ComputeGlobalResultsMsg;
  char *b = (char*)ptr;

  int size, gsize;
  memcpy(&size, b, sizeof(int)); b += sizeof(int);
  m->aid.resize(size);
  memcpy(m->aid.begin(), b, size*sizeof(AtomID)); b += size*sizeof(AtomID);
  m->f.resize(size);
  memcpy(m->f.begin(), b, size*sizeof(Force)); b += size*sizeof(Force);
  memcpy(&gsize, b, sizeof(int)); b += sizeof(int);
  m->gforce.resize(gsize);
  memcpy(m->gforce.begin(), b, gsize*sizeof(Force)); b += gsize*sizeof(Force);
  memcpy(&(m->reconfig), b, sizeof(int)); b += sizeof(int);
  if ( m->reconfig ) {
    memcpy(&size, b, sizeof(int)); b += sizeof(int);
    m->newaid.resize(size);
    memcpy(m->newaid.begin(), b, size*sizeof(AtomID)); b += size*sizeof(AtomID);
    memcpy(&gsize, b, sizeof(int)); b += sizeof(int);
    m->newgdef.resize(gsize);
    memcpy(m->newgdef.begin(), b, gsize*sizeof(AtomID)); b += gsize*sizeof(AtomID);
  }
  CkFreeMsg(ptr);
  return m;
}



/***************************************************************************
 * RCS INFORMATION:
 *
 *	$RCSfile: ComputeGlobalMsgs.C,v $
 *	$Author: brunner $	$Locker:  $		$State: Exp $
 *	$Revision: 1.4 $	$Date: 1999/05/11 23:56:22 $
 *
 * REVISION HISTORY:
 *
 * $Log: ComputeGlobalMsgs.C,v $
 * Revision 1.4  1999/05/11 23:56:22  brunner
 * Changes for new charm version
 *
 * Revision 1.3  1998/10/24 19:57:24  jim
 * Eliminated warnings generated by g++ -Wall.
 *
 * Revision 1.2  1998/02/16 00:24:38  jim
 * Added atom group centers of mass to Tcl interface.
 *
 * Revision 1.1  1997/12/19 23:48:47  jim
 * Added Tcl interface for calculating forces.
 *
 *
 ***************************************************************************/
