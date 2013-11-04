/**
***  Copyright (c) 1995, 1996, 1997, 1998, 1999, 2000 by
***  The Board of Trustees of the University of Illinois.
***  All rights reserved.
**/
#include <stdio.h>
#include "DataExchanger.h"
#include "ProcessorPrivate.h"
#include "common.h"

#if CMK_HAS_PARTITION
#ifdef CmiMyPartitionSize
extern "C" {
  void setDefaultPartitionParams() {
    // if(!CmiMyNodeGlobal()) printf("NAMD setDefaultPartitionParams called\n");
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

  void replica_bcast(char *buf, int count, int root) {
    if ( root != 0 ) NAMD_bug("replica_bcast root must be zero");
    int rank = CmiMyPartition();
    int size = CmiNumPartitions();
    int i;
    for ( i=1; i<size; i*=2 );
    for ( i/=2; i>0; i/=2 ) {
      if ( rank & (i - 1) ) continue;
      if ( rank & i ) {
        int src = rank - i;
        CkPrintf("rank %d recv from %d\n", rank, src);
        DataMessage *recvMsg = NULL;
        replica_recv(&recvMsg, src, CkMyPe());
        if ( recvMsg == NULL ) NAMD_bug("recvMsg == NULL in replica_bcast");
        if ( recvMsg->size != count ) NAMD_bug("size != count in replica_bcast");
        memcpy(buf, recvMsg->data, count);
        CmiFree(recvMsg);
      } else {
        int dst = rank + i;
        if ( dst < size ) {
          CkPrintf("rank %d send to %d\n", rank, dst);
          replica_send(buf, count, dst, CkMyPe());
        }
      }
    }
  }

  void replica_min_double(double *dat, int count) {
    int rank = CmiMyPartition();
    int size = CmiNumPartitions();
    for ( int i=1; i<size; i*=2 ) {
      if ( rank & i ) {
        int dst = rank - i;
        CkPrintf("rank %d send to %d\n", rank, dst);
        replica_send((char*)dat, count * sizeof(double), dst, CkMyPe());
      } else {
        int src = rank + i;
        if ( src < size ) {
          CkPrintf("rank %d recv from %d\n", rank, src);
          DataMessage *recvMsg = NULL;
          replica_recv(&recvMsg, src, CkMyPe());
          if ( recvMsg == NULL ) NAMD_bug("recvMsg == NULL in replica_bcast");
          if ( recvMsg->size != count * sizeof(double) ) NAMD_bug("size != count in replica_min_double");
          double *rdat = new double[count];
          memcpy(rdat, recvMsg->data, count * sizeof(double));
          CmiFree(recvMsg);
          for ( int j=0; j<count; ++j ) {
            if ( rdat[j] < dat[j] ) dat[j] = rdat[j];
          }
          delete [] rdat;
        }
      }
      if ( rank & (2 * i - 1) ) break;
    }
    replica_bcast((char*)dat, count * sizeof(double), 0);
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
