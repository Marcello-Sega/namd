/**
***  Copyright (c) 1995, 1996, 1997, 1998, 1999, 2000 by
***  The Board of Trustees of the University of Illinois.
***  All rights reserved.
**/
#include <stdio.h>
#include "DataExchanger.h"
#include "ProcessorPrivate.h"

#if CMK_HAS_PARTITION
#ifdef CmiMyPartitionSize
extern "C" {
  void setDefaultPartitionParams() {
    // if(!CmiMyNodeGlobal()) CmiPrintf("NAMD setDefaultPartitionParams called\n");
    CmiSetPartitionScheme(3);  // recursive bisection
  }
}
#endif
#endif

CpvDeclare(int, breakScheduler);

//functions to receive and invoke chare's entry methods
extern "C" {
  void packSend(int dst, int dstPart, char *data, int size, int handler) {
    int msgsize = sizeof(DataMessage) + size;
    DataMessage *dmsg = (DataMessage *)CmiAlloc(msgsize);
    dmsg->setMessage(data,CkMyPe(),CmiMyPartition(),size,handler);
#if CMK_HAS_PARTITION
    CmiInterSyncSendAndFree(dst,dstPart,msgsize,(char*)dmsg);
#else
    CmiSyncSendAndFree(dst,msgsize,(char*)dmsg);
#endif
  }

  void recvData(DataMessage *dmsg) {
    Pointer msg(dmsg);
    CPROXY_DE(CkpvAccess(BOCclass_group).dataExchanger)[CkMyPe()].recv_data(msg);
  }

  void recvAck(DataMessage *dmsg) {
    CPROXY_DE(CkpvAccess(BOCclass_group).dataExchanger)[CkMyPe()].recv_ack();
    CmiFree(dmsg);
  }

  void recvBcast(DataMessage *dmsg) {
    CPROXY_DE(CkpvAccess(BOCclass_group).dataExchanger)[CkMyPe()].recv_bcast();
    CmiFree(dmsg);
  }

  void recvRed(DataMessage *dmsg) {
    CPROXY_DE(CkpvAccess(BOCclass_group).dataExchanger)[CkMyPe()].recv_red();
    CmiFree(dmsg);
  }

  void replica_send(char *sndbuf, int sendcount, int destPart, int destPE) {
    Pointer sendPointer(sndbuf);
    CPROXY_DE(CkpvAccess(BOCclass_group).dataExchanger)[CkMyPe()].send(sendPointer,sendcount,destPart,destPE); 
    CpvAccess(breakScheduler) = 0;
    while(!CpvAccess(breakScheduler)) CsdSchedulePoll();
  }

  void replica_recv(DataMessage **precvMsg, int srcPart, int srcPE) {
    Pointer recvPointer((char *) precvMsg);
    CPROXY_DE(CkpvAccess(BOCclass_group).dataExchanger)[CkMyPe()].recv(recvPointer,srcPart,srcPE);
    CpvAccess(breakScheduler) = 0;
    while(!CpvAccess(breakScheduler)) CsdSchedulePoll();
  }

  void replica_sendRecv(char *sndbuf, int sendcount, int destPart, int destPE, DataMessage **precvMsg, int srcPart, int srcPE)  {
    Pointer sendPointer(sndbuf);
    Pointer recvPointer((char *) precvMsg);
    CPROXY_DE(CkpvAccess(BOCclass_group).dataExchanger)[CkMyPe()].sendRecv(sendPointer,sendcount,destPart,destPE,recvPointer,srcPart,srcPE);
    CpvAccess(breakScheduler) = 0;
    while(!CpvAccess(breakScheduler)) CsdSchedulePoll();
  }

  void replica_barrier() {
    CPROXY_DE(CkpvAccess(BOCclass_group).dataExchanger)[CkMyPe()].barrier();
    CpvAccess(breakScheduler) = 0;
    while(!CpvAccess(breakScheduler)) CsdSchedulePoll();
  }
} //endof extern C

//======================================================================
// Public functions
//----------------------------------------------------------------------
DataExchanger::DataExchanger()
{
  CpvInitialize(int, breakScheduler);
  CpvAccess(breakScheduler) = 1;
  if(CmiMyPartition() == 0) 
    parent = -1;
  else 
    parent = (CmiMyPartition()+1)/TREE_WIDTH - 1;
  firstChild = (CmiMyPartition()+1)*TREE_WIDTH - 1;
  numChildren = CmiNumPartitions() - firstChild;
  if(numChildren > TREE_WIDTH)
    numChildren = TREE_WIDTH;
  recv_data_idx = CmiRegisterHandler((CmiHandler)recvData);
  recv_ack_idx = CmiRegisterHandler((CmiHandler)recvAck);
  recv_red_idx = CmiRegisterHandler((CmiHandler)recvRed);
  recv_bcast_idx = CmiRegisterHandler((CmiHandler)recvBcast);
}

//----------------------------------------------------------------------
DataExchanger::~DataExchanger(void)
{ }


#include "DataExchanger.def.h"
