/**
***  Copyright (c) 1995, 1996, 1997, 1998, 1999, 2000 by
***  The Board of Trustees of the University of Illinois.
***  All rights reserved.
**/

#include "ComputeGlobalMsgs.h"
#include "packmsg.h"

//#define DEBUGM
#define MIN_DEBUG_LEVEL 3
#include "Debug.h"

#include "ComputeMgr.decl.h"


// CONFIG MESSAGE

ComputeGlobalConfigMsg::ComputeGlobalConfigMsg(void) : tag(-1) { 
}

ComputeGlobalConfigMsg::~ComputeGlobalConfigMsg(void) { 
}

PACK_MSG(ComputeGlobalConfigMsg,
  PACK(tag);
  PACK_RESIZE(aid);
  PACK_RESIZE(gdef);
)


// DATA MESSAGE

ComputeGlobalDataMsg::ComputeGlobalDataMsg(void) : tag(-1) { 
}

ComputeGlobalDataMsg::~ComputeGlobalDataMsg(void) { 
}

PACK_MSG(ComputeGlobalDataMsg,
  PACK(tag);
  PACK_RESIZE(aid);
  PACK_RESIZE(p);
  PACK_RESIZE(gcom);
)


// RESULTS MESSAGE

ComputeGlobalResultsMsg::ComputeGlobalResultsMsg(void) : tag(-1) { 
  reconfig = 0;
}

ComputeGlobalResultsMsg::~ComputeGlobalResultsMsg(void) { 
}

PACK_MSG(ComputeGlobalResultsMsg,
  PACK(tag);
  PACK_RESIZE(aid);
  PACK_RESIZE(f);
  PACK_RESIZE(gforce);
  PACK(reconfig);
  if ( packmsg_msg->reconfig ) {
    PACK_RESIZE(newaid);
    PACK_RESIZE(newgdef);
  }
)

