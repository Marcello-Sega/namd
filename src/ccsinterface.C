/**
***  Copyright (c) 1995, 1996, 1997, 1998, 1999, 2000 by
***  The Board of Trustees of the University of Illinois.
***  All rights reserved.
**/

#include "ccsinterface.h"

#ifdef NAMDCCS
#include <stdlib.h>

extern unsigned int appletIP;
extern unsigned int appletPort;
CpvDeclare(int, CApplicationDataCollectionHandlerIndex);
int applicationCountMsgs;
char **applicationValueArray;

static void sendDataFunction(void)
{
  char *reply;
  int len = 0, i;

  for(i=0; i<CmiNumPes(); i++){
    len += (strlen((char*)(applicationValueArray[i]+
			   CmiMsgHeaderSizeBytes+sizeof(int)))+1);
    /* for the spaces in between */
  }
  len+=8; /* for 'robert ' and the \0 at the end */

  reply = (char *)malloc(len * sizeof(char));
  strcpy(reply, "robert ");

  for(i=0; i<CmiNumPes(); i++){
    strcat(reply, (applicationValueArray[i] + 
		   CmiMsgHeaderSizeBytes + sizeof(int)));
    strcat(reply, " ");
  }

  /* Do the CcsSendReply */
#if CMK_USE_PERSISTENT_CCS
  CcsSendReplyFd(appletIP, appletPort, strlen(reply) + 1, reply);
#else
  CcsSendReply(appletIP, appletPort, strlen(reply) + 1, reply);
#endif
  CmiPrintf("Reply = %s\n", reply);
  free(reply);

  /* Free applicationValueArray contents */
  for(i = 0; i < CmiNumPes(); i++){
    CmiFree(applicationValueArray[i]);
    applicationValueArray[i] = 0;
  }

  applicationCountMsgs = 0;
}

void CApplicationDataCollectionHandler(char *msg){
  int src;
  int value;
  char *prev;

  if(CmiMyPe() != 0){
    CmiAbort("Wrong processor....\n");
  }
  src = ((int *)(msg + CmiMsgHeaderSizeBytes))[0];
  CmiGrabBuffer((void **)&msg);
  prev = applicationValueArray[src]; /* Previous value, ideally 0 */
  applicationValueArray[src] = (msg);
  if(prev == 0) applicationCountMsgs++;
  else CmiFree(prev);

  if(applicationCountMsgs == CmiNumPes()){
    sendDataFunction();
  }
}

extern "C" void CApplicationDepositData(char *data)
{
  char *msg;
  int msgSize;
  int i;

  if(appletIP == 0) {
    return; 
  }

  msgSize = (strlen(data)+1) + sizeof(int) + CmiMsgHeaderSizeBytes;
  msg = (char *)CmiAlloc(msgSize);
  ((int *)(msg + CmiMsgHeaderSizeBytes))[0] = CmiMyPe();
  strcpy(msg + CmiMsgHeaderSizeBytes + sizeof(int), data);
  CmiSetHandler(msg, CpvAccess(CApplicationDataCollectionHandlerIndex));
  CmiSyncSendAndFree(0, msgSize, msg);
}

extern "C" void CApplicationDepositNode0Data(char *data)
{
  char *reply;
  int len;

  if(appletIP == 0) {
    return; 
  }
  
  len = strlen(data) + 8; /* for the 'robert ' and '\0' */
  reply = (char *)malloc(len * sizeof(char));
  strcpy(reply, "namdpr ");
  strcat(reply, data);
  
  /* Do the CcsSendReply */
#if CMK_USE_PERSISTENT_CCS
  CcsSendReplyFd(appletIP, appletPort, strlen(reply) + 1, reply);
#else
  CcsSendReply(appletIP, appletPort, strlen(reply) + 1, reply);
#endif
  CmiPrintf("Reply = %s\n", reply);
  free(reply);
}

void CApplicationInit(void)
{
  int i;

  CpvInitialize(int, CApplicationDataCollectionHandlerIndex);
  CpvAccess(CApplicationDataCollectionHandlerIndex) =
    CmiRegisterHandler((CmiHandler)CApplicationDataCollectionHandler);

  applicationValueArray = (char **)malloc(sizeof(char *) * CmiNumPes());
  for(i = 0; i < CmiNumPes(); i++)
    applicationValueArray[i] = 0;
}

#endif
