/**
 ***  Copyright (c) 1995, 1996, 1997, 1998, 1999, 2000 by
 ***  The Board of Trustees of the University of Illinois.
 ***  All rights reserved.
 **/

#ifndef _DATAEXCHANGER_H
#define _DATAEXCHANGER_H

#include "charm++.h"
#include "main.h"
#include "DataExchanger.decl.h"

#define CPROXY_DE(x) ((CProxy_DataExchanger)(x))
CpvExtern(int, breakScheduler);

class DataMessage {
  public:
  char core[CmiMsgHeaderSizeBytes];
  int src, srcPart;
  int size;
  char data[1];

  void setMessage(char *_data, int _src, int _srcPart, int _size, int _handler) {
    src = _src; srcPart = _srcPart;
    size = _size;
    memcpy(data,_data,size);
    CmiSetHandler(core,_handler);
  }
};

class DataExchanger : public CBase_DataExchanger
{
  public:
    DataExchanger_SDAG_CODE;
    DataExchanger();
    ~DataExchanger(void);
  
    int loop, recvred, sendbcast;
    enum{ TREE_WIDTH=2};
    int numChildren, firstChild, parent;

    //message handlers
    int recv_data_idx;
    int recv_ack_idx;
    int recv_bcast_idx;
    int recv_red_idx;
};

extern "C" {
void packSend(int dest, int partition, char *data, int size, int handler);
void recvData(DataMessage *dmsg); 
void recvAck(DataMessage *dmsg); 

void replica_send(char *sndbuf, int sendcount, int destPart, int destPE);
void replica_sendRecv(char *sndbuf, int sendcount, int destPart, int destPE, DataMessage **precvMsg, int srcPart, int srcPE);
void replica_recv(DataMessage **precvMsg, int srcPart, int srcPE);
void replica_barrier();
}
#endif
