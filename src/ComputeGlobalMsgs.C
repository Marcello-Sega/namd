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
#include "packmsg.h"

//#define DEBUGM
#define MIN_DEBUG_LEVEL 3
#include "Debug.h"

#include "ComputeMgr.decl.h"


// CONFIG MESSAGE

ComputeGlobalConfigMsg::ComputeGlobalConfigMsg(void) { 
}

ComputeGlobalConfigMsg::~ComputeGlobalConfigMsg(void) { 
}

PACK_MSG(ComputeGlobalConfigMsg,
  PACK_RESIZE(aid);
  PACK_RESIZE(gdef);
)


// DATA MESSAGE

ComputeGlobalDataMsg::ComputeGlobalDataMsg(void) { 
}

ComputeGlobalDataMsg::~ComputeGlobalDataMsg(void) { 
}

PACK_MSG(ComputeGlobalDataMsg,
  PACK_RESIZE(aid);
  PACK_RESIZE(p);
  PACK_RESIZE(gcom);
)


// RESULTS MESSAGE

ComputeGlobalResultsMsg::ComputeGlobalResultsMsg(void) { 
  reconfig = 0;
}

ComputeGlobalResultsMsg::~ComputeGlobalResultsMsg(void) { 
}

PACK_MSG(ComputeGlobalResultsMsg,
  PACK_RESIZE(aid);
  PACK_RESIZE(f);
  PACK_RESIZE(gforce);
  PACK(reconfig);
  if ( packmsg_msg->reconfig ) {
    PACK_RESIZE(newaid);
    PACK_RESIZE(newgdef);
  }
)

