/**
***  Copyright (c) 1995, 1996, 1997, 1998, 1999, 2000 by
***  The Board of Trustees of the University of Illinois.
***  All rights reserved.
**/

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
  energy = 0.;
  for ( int i=0; i<6; ++i ) { virial[i] = 0.; }
  start = 0;
  q_len = 0;
  q_arr = 0;
  will_delete_array = 1;
}

ComputePmeResultsMsg::~ComputePmeResultsMsg(void) { 
  if ( will_delete_array ) { delete [] q_arr; }
}

PACK_MSG(ComputePmeResultsMsg,
  PACK(energy);
  PACK(virial);
  PACK(start);
  PACK_AND_NEW_ARRAY(q_arr,q_len);
)

