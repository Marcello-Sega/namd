/**
***  Copyright (c) 1995, 1996, 1997, 1998, 1999, 2000 by
***  The Board of Trustees of the University of Illinois.
***  All rights reserved.
**/

#include "SMDMsgs.h"

// #define DEBUGM
#define MIN_DEBUG_LEVEL 3
#include "Debug.h"

// DATA MESSAGE

void * SMDDataMsg::pack (SMDDataMsg *m) {
  int length = 2 * sizeof(int) + 4 * sizeof(Vector);

  char *buffer;
  char *b = buffer = (char*)CkAllocBuffer(m,length);

  memcpy(b, &(m->curTime), sizeof(int)); b += sizeof(int);
  memcpy(b, &(m->timeStamp), sizeof(int)); b += sizeof(int);
  memcpy(b, &(m->direction), sizeof(Vector)); b += sizeof(Vector);
  memcpy(b, &(m->refPos), sizeof(Vector)); b += sizeof(Vector);
  memcpy(b, &(m->atomPosVmin), sizeof(Vector)); b += sizeof(Vector);
  memcpy(b, &(m->atomPosVmax), sizeof(Vector)); b += sizeof(Vector);

  delete m;
  return buffer;
}

SMDDataMsg* SMDDataMsg::unpack (void *ptr) {
  void *_ptr = CkAllocBuffer(ptr, sizeof(SMDDataMsg));
  SMDDataMsg* m = new (_ptr) SMDDataMsg;
  char *b = (char*)ptr;

  memcpy(&(m->curTime), b, sizeof(int)); b += sizeof(int);
  memcpy(&(m->timeStamp), b, sizeof(int)); b += sizeof(int);
  memcpy(&(m->direction), b, sizeof(Vector)); b += sizeof(Vector);
  memcpy(&(m->refPos), b, sizeof(Vector)); b += sizeof(Vector);
  memcpy(&(m->atomPosVmin), b, sizeof(Vector)); b += sizeof(Vector);
  memcpy(&(m->atomPosVmax), b, sizeof(Vector)); b += sizeof(Vector);

  CkFreeMsg(ptr);
  return m;
}

