/***************************************************************************/
/*         (C) Copyright 1996,1997 The Board of Trustees of the            */
/*                          University of Illinois                         */
/*                           All Rights Reserved                           */
/***************************************************************************/
/***************************************************************************
 * DESCRIPTION:  Messages needed for ComputeDPME operation.
 *		
 ***************************************************************************/

static char ident[] = "@(#)$Header: /home/cvs/namd/cvsroot/namd2/src/ComputeDPMEMsgs.C,v 1.1 1998/04/10 04:15:58 jim Exp $";

#include "ComputeDPMEMsgs.h"

//#define DEBUGM
#define MIN_DEBUG_LEVEL 3
#include "Debug.h"

#include "ComputeMgr.top.h"

#ifdef DPME
#include "dpme2.h"
#else
#define Pme2Particle char;
#define PmeVector char;
#endif

// DATA MESSAGE

ComputeDPMEDataMsg::ComputeDPMEDataMsg(void) { 
  numParticles = 0;
  particles = 0;
}

ComputeDPMEDataMsg::~ComputeDPMEDataMsg(void) { 
  delete [] particles;
}

void * ComputeDPMEDataMsg::pack (int *length) {
  *length = 2 * sizeof(int) + numParticles * sizeof(Pme2Particle);

  char *buffer;
  char *b = buffer = (char*)new_packbuffer(this,*length);

  memcpy(b, &node, sizeof(int)); b += sizeof(int);
  memcpy(b, &numParticles, sizeof(int)); b += sizeof(int);
  memcpy(b, particles, numParticles*sizeof(Pme2Particle));
  b += numParticles*sizeof(Pme2Particle);

  this->~ComputeDPMEDataMsg();
  return buffer;
}

void ComputeDPMEDataMsg::unpack (void *in) {
  new((void*)this) ComputeDPMEDataMsg;
  char *b = (char*)in;

  memcpy(&node, b, sizeof(int)); b += sizeof(int);
  memcpy(&numParticles, b, sizeof(int)); b += sizeof(int);
  particles = new Pme2Particle[numParticles];
  memcpy(particles, b, numParticles*sizeof(Pme2Particle));
  b += numParticles*sizeof(Pme2Particle);

  // DO NOT delete void *in - this is done by Charm
}


// RESULTS MESSAGE

ComputeDPMEResultsMsg::ComputeDPMEResultsMsg(void) { 
  numParticles = 0;
  forces = 0;
}

ComputeDPMEResultsMsg::~ComputeDPMEResultsMsg(void) { 
  delete [] forces;
}

void * ComputeDPMEResultsMsg::pack (int *length) {
  *length = 2 * sizeof(int) + numParticles * sizeof(Pme2Particle);

  char *buffer;
  char *b = buffer = (char*)new_packbuffer(this,*length);

  memcpy(b, &node, sizeof(int)); b += sizeof(int);
  memcpy(b, &numParticles, sizeof(int)); b += sizeof(int);
  memcpy(b, forces, numParticles*sizeof(PmeVector));
  b += numParticles*sizeof(PmeVector);

  this->~ComputeDPMEResultsMsg();
  return buffer;
}

void ComputeDPMEResultsMsg::unpack (void *in) {
  new((void*)this) ComputeDPMEResultsMsg;
  char *b = (char*)in;

  memcpy(&node, b, sizeof(int)); b += sizeof(int);
  memcpy(&numParticles, b, sizeof(int)); b += sizeof(int);
  forces = new PmeVector[numParticles];
  memcpy(forces, b, numParticles*sizeof(PmeVector));
  b += numParticles*sizeof(PmeVector);

  // DO NOT delete void *in - this is done by Charm
}



/***************************************************************************
 * RCS INFORMATION:
 *
 *	$RCSfile: ComputeDPMEMsgs.C,v $
 *	$Author: jim $	$Locker:  $		$State: Exp $
 *	$Revision: 1.1 $	$Date: 1998/04/10 04:15:58 $
 *
 * REVISION HISTORY:
 *
 * $Log: ComputeDPMEMsgs.C,v $
 * Revision 1.1  1998/04/10 04:15:58  jim
 * Finished incorporating DPME.
 *
 *
 ***************************************************************************/
