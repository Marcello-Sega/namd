/**
***  Copyright (c) 1995, 1996, 1997, 1998, 1999, 2000 by
***  The Board of Trustees of the University of Illinois.
***  All rights reserved.
**/

#include <string.h>
#include <stdlib.h>
#include "Communicate.h"
#include "MStream.h"
#include "converse.h"

CpvStaticDeclare(CmmTable, CsmMessages);

static void CsmHandler(void *msg)
{
  if ( CmiMyRank() ) NAMD_bug("CsmHandler called on non-rank-zero Pe\n");
  // if ( CmiMyRank() ) { CmiFree(msg); return; }
  // get start of user message
  int *m = (int *) ((char *)msg+CmiMsgHeaderSizeBytes);
  // sending node  & tag act as tags
  CmmPut(CpvAccess(CsmMessages), 2, m, msg);
}

CmiGroup rankZeroPes;

Communicate::Communicate(void) 
{
  CpvInitialize(CmmTable, CsmMessages);
  CsmHandlerIndex = CmiRegisterHandler((CmiHandler) CsmHandler);
  CpvAccess(CsmMessages) = CmmNew();

  if ( CmiMyPe() == 0 ) {
    int *pes = new int[CmiNumPes()];
    int npes = 0;
    for ( int i = 1; i < CmiNumPes(); ++i ) {
      if ( CmiRankOf(i) == 0 ) { pes[npes++] = i; }
    }
    for ( int j = 0; j < npes; ++j ) { CmiPrintf("%d\n",pes[j]); }
    rankZeroPes = CmiEstablishGroup(npes,pes);
    delete [] pes;
  }
}


Communicate::~Communicate(void) 
{
  // do nothing
}

MIStream *Communicate::newInputStream(int PE, int tag)
{
  MIStream *st = new MIStream(this, PE, tag);
  return st;
}

MOStream *Communicate::newOutputStream(int PE, int tag, unsigned int bufSize)
{
  MOStream *st = new MOStream(this, PE, tag, bufSize);
  return st;
}

void *Communicate::getMessage(int PE, int tag)
{
  if ( CmiMyRank() ) NAMD_bug("Communicate::getMessage called on non-rank-zero Pe\n");
  int itag[2], rtag[2];
  void *msg;

  itag[0] = (PE==(-1)) ? (CmmWildCard) : PE;
  itag[1] = (tag==(-1)) ? (CmmWildCard) : tag;
  while((msg=CmmGet(CpvAccess(CsmMessages),2,itag,rtag))==0) {
    CmiDeliverMsgs(0);
  }
  return msg;
}

void Communicate::sendMessage(int PE, void *msg, int size)
{
  if ( CmiMyPe() ) NAMD_bug("Communicate::sendMessage not from Pe 0");
  CmiSetHandler(msg, CsmHandlerIndex);
  switch(PE) {
    case ALL:
      NAMD_bug("Unexpected Communicate::sendMessage(ALL,...)");
      //CmiSyncBroadcastAll(size, (char *)msg);
      break;
    case ALLBUTME:
      CmiSyncMulticast(rankZeroPes, size, (char *)msg);
      // CmiSyncBroadcast(size, (char *)msg);
      break;
    default:
      NAMD_bug("Unexpected Communicate::sendMessage(PEL,...)");
      //CmiSyncSend(PE, size, (char *)msg);
      break;
  }
}
