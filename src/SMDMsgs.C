/***************************************************************************/
/*         (C) Copyright 1996,1997 The Board of Trustees of the            */
/*                          University of Illinois                         */
/*                           All Rights Reserved                           */
/***************************************************************************/
/***************************************************************************
 * DESCRIPTION:  Messages needed for sending SMD data.
 *		
 ***************************************************************************/

#include "SMDMsgs.h"

// #define DEBUGM
#define MIN_DEBUG_LEVEL 3
#include "Debug.h"

// DATA MESSAGE

void * SMDDataMsg::pack (int *length) {
  *length = 2 * sizeof(int) + 4 * sizeof(Vector);

  char *buffer;
  char *b = buffer = (char*)new_packbuffer(this,*length);

  memcpy(b, &curTime, sizeof(int)); b += sizeof(int);
  memcpy(b, &timeStamp, sizeof(int)); b += sizeof(int);
  memcpy(b, &direction, sizeof(Vector)); b += sizeof(Vector);
  memcpy(b, &refPos, sizeof(Vector)); b += sizeof(Vector);
  memcpy(b, &atomPosVmin, sizeof(Vector)); b += sizeof(Vector);
  memcpy(b, &atomPosVmax, sizeof(Vector)); b += sizeof(Vector);

  this->~SMDDataMsg();
  return buffer;
}

void SMDDataMsg::unpack (void *in) {
  new((void*)this) SMDDataMsg;
  char *b = (char*)in;

  memcpy(&curTime, b, sizeof(int)); b += sizeof(int);
  memcpy(&timeStamp, b, sizeof(int)); b += sizeof(int);
  memcpy(&direction, b, sizeof(Vector)); b += sizeof(Vector);
  memcpy(&refPos, b, sizeof(Vector)); b += sizeof(Vector);
  memcpy(&atomPosVmin, b, sizeof(Vector)); b += sizeof(Vector);
  memcpy(&atomPosVmax, b, sizeof(Vector)); b += sizeof(Vector);

  // DO NOT delete void *in - this is done by Charm
}

