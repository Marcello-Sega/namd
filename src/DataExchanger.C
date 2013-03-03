/**
***  Copyright (c) 1995, 1996, 1997, 1998, 1999, 2000 by
***  The Board of Trustees of the University of Illinois.
***  All rights reserved.
**/
#include <stdio.h>
#include "DataExchanger.h"
#include "ProcessorPrivate.h"

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
    CsdScheduler(-1);
  }

  void replica_recv(char *precvMsg, int srcPart, int srcPE) {
    Pointer recvPointer(precvMsg);
    CPROXY_DE(CkpvAccess(BOCclass_group).dataExchanger)[CkMyPe()].recv(recvPointer,srcPart,srcPE);
    CsdScheduler(-1);
  }

  void replica_sendRecv(char *sndbuf, int sendcount, int destPart, int destPE, char *precvMsg, int srcPart, int srcPE)  {
    Pointer sendPointer(sndbuf);
    Pointer recvPointer(precvMsg);
    CPROXY_DE(CkpvAccess(BOCclass_group).dataExchanger)[CkMyPe()].sendRecv(sendPointer,sendcount,destPart,destPE,recvPointer,srcPart,srcPE);
    CsdScheduler(-1);
  }

  void replica_barrier() {
    CPROXY_DE(CkpvAccess(BOCclass_group).dataExchanger)[CkMyPe()].barrier();
    CsdScheduler(-1);
  }
} //endof extern C

//======================================================================
// Public functions
//----------------------------------------------------------------------
DataExchanger::DataExchanger()
{
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
