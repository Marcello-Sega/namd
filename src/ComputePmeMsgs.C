/***************************************************************************/
/*         (C) Copyright 1996,1997 The Board of Trustees of the            */
/*                          University of Illinois                         */
/*                           All Rights Reserved                           */
/***************************************************************************/
/***************************************************************************
 * DESCRIPTION:  Messages needed for ComputePme operation.
 *		
 ***************************************************************************/

#include "ComputePmeMsgs.h"
#include "packmsg.h"

//#define DEBUGM
#define MIN_DEBUG_LEVEL 3
#include "Debug.h"

#include "ComputeMgr.decl.h"

#include "PmeBase.h"


// DATA MESSAGE

ComputePmeDataMsg::ComputePmeDataMsg(void) { 
  numParticles = 0;
  particles = 0;
}

ComputePmeDataMsg::~ComputePmeDataMsg(void) { 
  delete [] particles;
}

PACK_MSG(ComputePmeDataMsg,
  PACK(node);
  PACK_AND_NEW_ARRAY(particles,numParticles);
)


// RESULTS MESSAGE

ComputePmeResultsMsg::ComputePmeResultsMsg(void) { 
  numParticles = 0;
  forces = 0;
}

ComputePmeResultsMsg::~ComputePmeResultsMsg(void) { 
  delete [] forces;
}

PACK_MSG(ComputePmeResultsMsg,
  PACK(node);
  PACK_AND_NEW_ARRAY(forces,numParticles);
)

