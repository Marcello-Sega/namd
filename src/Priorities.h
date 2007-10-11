/**
***  Copyright (c) 1995, 1996, 1997, 1998, 1999, 2000 by
***  The Board of Trustees of the University of Illinois.
***  All rights reserved.
**/

#ifndef PRIORITIES_H
#define PRIORITIES_H

// pass PRIORITY_SIZE as last argument to new when allocating message
// e.g.,  MyMsg *msg = new (len1, len2, PRIORITY_SIZE) MyMsg;

#define PRIORITY_SIZE ((int) sizeof(int)*8)

// total priority is sequence + type + patch
// always use lowest (most urgent) patch priority applicable

#define SET_PRIORITY(MSG,SEQ,PRIO) { \
  CkSetQueueing(MSG, CK_QUEUEING_IFIFO); \
  *((int*) CkPriorityPtr(MSG)) = (((SEQ)&0xffff)<<15) + (PRIO); }
// sequence priority is 16 bits shifted by 15 to leave sign bit 0

// patch priorities in range [1,252]
// reserve 0 for tuples, use prime to randomize neighbors
#define PATCH_PRIORITY(PID) (((PID)%251)+1)

// the following are in order of decreasing urgency

#define PME_GRID_PRIORITY (1<<8)
#define PME_TRANS_PRIORITY (PME_GRID_PRIORITY+1)
#define PME_TRANS2_PRIORITY (PME_GRID_PRIORITY+2)
#define PME_UNTRANS_PRIORITY (PME_GRID_PRIORITY+3)
#define PME_UNTRANS2_PRIORITY (PME_GRID_PRIORITY+4)

#define PROXY_DATA_PRIORITY (2<<8)

#define COMPUTE_PROXY_PRIORITY (3<<8)

#define PROXY_RESULTS_PRIORITY (4<<8)

// pme ungrid messages only affect home patches
#define PME_UNGRID_PRIORITY (5<<8)

#define COMPUTE_HOME_PRIORITY (6<<8)

#endif // PRIORITIES_H

